#define	JEMALLOC_TCACHE_C_
#include "jemalloc/internal/jemalloc_internal.h"

/******************************************************************************/
/* Data. */

bool	opt_tcache = true;
ssize_t	opt_lg_tcache_max = LG_TCACHE_MAXCLASS_DEFAULT;

tcache_bin_info_t	*tcache_bin_info;
extern arena_cache_bin_info_t	*arena_cache_bin_info;

unsigned		tcache_stack_nelms; /* Total stack elms per tcache. */

unsigned		nhbins;
size_t			tcache_maxclass;

tcaches_t		*tcaches;

/* Index of first element within tcaches that has never been used. */
static unsigned		tcaches_past;

/* Head of singly linked list tracking available tcaches elements. */
static tcaches_t	*tcaches_avail;

/******************************************************************************/

size_t
tcache_salloc(tsdn_t *tsdn, const void *ptr)
{

	return (arena_salloc(tsdn, ptr, false));
}

void
tcache_event_hard(tsd_t *tsd, tcache_t *tcache)
{
	szind_t binind = tcache->next_gc_bin;
	tcache_bin_t *tbin = &tcache->tbins[binind];
	tcache_bin_info_t *tbin_info = &tcache_bin_info[binind];

	if (tbin->low_water > 0) {
		/*
		 * Flush (ceiling) 3/4 of the objects below the low water mark.
		 */
		unsigned flush_rem = tbin->ncached - tbin->low_water
		    + (tbin->low_water >> 2);
		tcache_bin_flush(tsd, tcache, tbin, binind, flush_rem,
		    binind < NBINS ? false : true);

		/*
		 * Reduce fill count by 2X.  Limit lg_fill_div such that the
		 * fill count is always at least 1.
		 */
		if ((tbin_info->ncached_max >> (tbin->lg_fill_div+1)) >= 1)
			tbin->lg_fill_div++;
	} else if (tbin->low_water < 0) {
		/*
		 * Increase fill count by 2X.  Make sure lg_fill_div stays
		 * greater than 0.
		 */
		if (tbin->lg_fill_div > 1)
			tbin->lg_fill_div--;
	}
	tbin->low_water = tbin->ncached;

	tcache->next_gc_bin++;
	if (tcache->next_gc_bin == nhbins)
		tcache->next_gc_bin = 0;

	if (config_acache && opt_acache) {
		arena_cache_gc(tsd_tsdn(tsd), arena_choose(tsd, NULL));
	}
}

void *
tcache_alloc_small_hard(tsdn_t *tsdn, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, szind_t binind, bool *tcache_success)
{
	void *ret;

	if (!config_acache || !opt_acache ||
	    !arena_cache_alloc_small(tsdn, arena, tcache, tbin, binind)) {
		/* Acache not available. Fill from arena directly. */
		arena_tcache_fill_small(tsdn, arena, tbin, binind, config_prof ?
		    tcache->prof_accumbytes : 0);
		if (config_prof)
			tcache->prof_accumbytes = 0;
	}

	ret = tcache_alloc_easy(tbin, tcache_success);

	return (ret);
}

void
tcache_bin_flush(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    szind_t binind, unsigned rem, const bool is_large)
{
	arena_t *arena;
	unsigned nflush;

	if (is_large)
		assert(binind < nhbins);
	else
		assert(binind < NBINS);
	assert(tbin->ncached && rem <= tbin->ncached);

	arena = arena_choose(tsd, NULL);
	assert(arena != NULL);

	nflush = tbin->ncached - rem;
	arena_cache_dalloc(tsd_tsdn(tsd), arena, tbin->avail - nflush, nflush,
	    binind, tbin->tstats.nrequests, is_large);

	if (config_prof) {
		if (arena_prof_accum(tsd_tsdn(tsd), arena, tcache->prof_accumbytes))
			prof_idump(tsd_tsdn(tsd));
		tcache->prof_accumbytes = 0;
	}
	tbin->tstats.nrequests = 0;

	memmove(tbin->avail - rem, tbin->avail - tbin->ncached, rem *
	    sizeof(void *));
	tbin->ncached = rem;
	if ((int)tbin->ncached < tbin->low_water)
		tbin->low_water = tbin->ncached;
}

static void
tcache_arena_associate(tsdn_t *tsdn, tcache_t *tcache, arena_t *arena)
{

	if (config_stats) {
		/* Link into list of extant tcaches. */
		malloc_mutex_lock(tsdn, &arena->lock);
		ql_elm_new(tcache, link);
		ql_tail_insert(&arena->tcache_ql, tcache, link);
		malloc_mutex_unlock(tsdn, &arena->lock);
	}
}

static void
tcache_arena_dissociate(tsdn_t *tsdn, tcache_t *tcache, arena_t *arena)
{

	if (config_stats) {
		/* Unlink from list of extant tcaches. */
		malloc_mutex_lock(tsdn, &arena->lock);
		if (config_debug) {
			bool in_ql = false;
			tcache_t *iter;
			ql_foreach(iter, &arena->tcache_ql, link) {
				if (iter == tcache) {
					in_ql = true;
					break;
				}
			}
			assert(in_ql);
		}
		ql_remove(&arena->tcache_ql, tcache, link);
		tcache_stats_merge(tsdn, tcache, arena);
		malloc_mutex_unlock(tsdn, &arena->lock);
	}
}

void
tcache_arena_reassociate(tsdn_t *tsdn, tcache_t *tcache, arena_t *oldarena,
    arena_t *newarena)
{

	tcache_arena_dissociate(tsdn, tcache, oldarena);
	tcache_arena_associate(tsdn, tcache, newarena);
}

tcache_t *
tcache_get_hard(tsd_t *tsd)
{
	arena_t *arena;

	if (!tcache_enabled_get()) {
		if (tsd_nominal(tsd))
			tcache_enabled_set(false); /* Memoize. */
		return (NULL);
	}
	arena = arena_choose(tsd, NULL);
	if (unlikely(arena == NULL))
		return (NULL);

	return (tcache_create(tsd_tsdn(tsd), arena));
}

tcache_t *
tcache_create(tsdn_t *tsdn, arena_t *arena)
{
	tcache_t *tcache;
	size_t size, stack_offset;
	unsigned i;

	size = offsetof(tcache_t, tbins) + (sizeof(tcache_bin_t) * nhbins);
	/* Naturally align the pointer stacks. */
	size = PTR_CEILING(size);
	stack_offset = size;
	size += tcache_stack_nelms * sizeof(void *);
	/* Avoid false cacheline sharing. */
	size = sa2u(size, CACHELINE);

	tcache = ipallocztm(tsdn, size, CACHELINE, true, NULL, true,
	    arena_get(TSDN_NULL, 0, true));
	if (tcache == NULL)
		return (NULL);

	tcache_arena_associate(tsdn, tcache, arena);

	ticker_init(&tcache->gc_ticker, TCACHE_GC_INCR);

	assert((TCACHE_NSLOTS_SMALL_MAX & 1U) == 0);
	for (i = 0; i < nhbins; i++) {
		tcache->tbins[i].lg_fill_div = 1;
		stack_offset += tcache_bin_info[i].ncached_max * sizeof(void *);
		/*
		 * avail points past the available space.  Allocations will
		 * access the slots toward higher addresses (for the benefit of
		 * prefetch).
		 */
		tcache->tbins[i].avail = (void **)((uintptr_t)tcache +
		    (uintptr_t)stack_offset);
	}

	return (tcache);
}

static void
tcache_destroy(tsd_t *tsd, tcache_t *tcache)
{
	arena_t *arena;
	unsigned i;

	arena = arena_choose(tsd, NULL);

	for (i = 0; i < NBINS; i++) {
		tcache_bin_t *tbin = &tcache->tbins[i];

		if (tbin->ncached)
			tcache_bin_flush(tsd, tcache, tbin, i, 0, false);

		if (config_stats && tbin->tstats.nrequests != 0) {
			arena_bin_t *bin = &arena->bins[i];
			malloc_mutex_lock(tsd_tsdn(tsd), &bin->lock);
			bin->stats.nrequests += tbin->tstats.nrequests;
			malloc_mutex_unlock(tsd_tsdn(tsd), &bin->lock);
		}
	}

	for (; i < nhbins; i++) {
		tcache_bin_t *tbin = &tcache->tbins[i];

		if (tbin->ncached)
			tcache_bin_flush(tsd, tcache, tbin, i, 0, true);

		if (config_stats && tbin->tstats.nrequests != 0) {
			malloc_mutex_lock(tsd_tsdn(tsd), &arena->lock);
			arena->stats.nrequests_large += tbin->tstats.nrequests;
			arena->stats.lstats[i - NBINS].nrequests +=
			    tbin->tstats.nrequests;
			malloc_mutex_unlock(tsd_tsdn(tsd), &arena->lock);
		}
	}

	if (opt_perCPU_arena) {
		/* Associated arena could have changed during flush. */
		arena = arena_choose(tsd, NULL);
	}
	tcache_arena_dissociate(tsd_tsdn(tsd), tcache, arena);

	if (config_prof && tcache->prof_accumbytes > 0 &&
	    arena_prof_accum(tsd_tsdn(tsd), arena, tcache->prof_accumbytes))
		prof_idump(tsd_tsdn(tsd));

	idalloctm(tsd_tsdn(tsd), tcache, NULL, true, true);
}

void
tcache_cleanup(tsd_t *tsd)
{
	tcache_t *tcache;

	if (!config_tcache)
		return;

	if ((tcache = tsd_tcache_get(tsd)) != NULL) {
		tcache_destroy(tsd, tcache);
		tsd_tcache_set(tsd, NULL);
	}
}

void
tcache_enabled_cleanup(tsd_t *tsd)
{

	/* Do nothing. */
}

void
tcache_stats_merge(tsdn_t *tsdn, tcache_t *tcache, arena_t *arena)
{
	unsigned i;

	cassert(config_stats);

	malloc_mutex_assert_owner(tsdn, &arena->lock);

	/* Merge and reset tcache stats. */
	for (i = 0; i < NBINS; i++) {
		arena_bin_t *bin = &arena->bins[i];
		tcache_bin_t *tbin = &tcache->tbins[i];
		malloc_mutex_lock(tsdn, &bin->lock);
		bin->stats.nrequests += tbin->tstats.nrequests;
		malloc_mutex_unlock(tsdn, &bin->lock);
		tbin->tstats.nrequests = 0;
	}

	for (; i < nhbins; i++) {
		malloc_large_stats_t *lstats = &arena->stats.lstats[i - NBINS];
		tcache_bin_t *tbin = &tcache->tbins[i];
		arena->stats.nrequests_large += tbin->tstats.nrequests;
		lstats->nrequests += tbin->tstats.nrequests;
		tbin->tstats.nrequests = 0;
	}
}

bool
tcaches_create(tsdn_t *tsdn, unsigned *r_ind)
{
	arena_t *arena;
	tcache_t *tcache;
	tcaches_t *elm;

	if (tcaches == NULL) {
		tcaches = base_alloc(tsdn, sizeof(tcache_t *) *
		    (MALLOCX_TCACHE_MAX+1));
		if (tcaches == NULL)
			return (true);
	}

	if (tcaches_avail == NULL && tcaches_past > MALLOCX_TCACHE_MAX)
		return (true);
	arena = arena_ichoose(tsdn, NULL);
	if (unlikely(arena == NULL))
		return (true);
	tcache = tcache_create(tsdn, arena);
	if (tcache == NULL)
		return (true);

	if (tcaches_avail != NULL) {
		elm = tcaches_avail;
		tcaches_avail = tcaches_avail->next;
		elm->tcache = tcache;
		*r_ind = (unsigned)(elm - tcaches);
	} else {
		elm = &tcaches[tcaches_past];
		elm->tcache = tcache;
		*r_ind = tcaches_past;
		tcaches_past++;
	}

	return (false);
}

static void
tcaches_elm_flush(tsd_t *tsd, tcaches_t *elm)
{

	if (elm->tcache == NULL)
		return;
	tcache_destroy(tsd, elm->tcache);
	elm->tcache = NULL;
}

void
tcaches_flush(tsd_t *tsd, unsigned ind)
{

	tcaches_elm_flush(tsd, &tcaches[ind]);
}

void
tcaches_destroy(tsd_t *tsd, unsigned ind)
{
	tcaches_t *elm = &tcaches[ind];
	tcaches_elm_flush(tsd, elm);
	elm->next = tcaches_avail;
	tcaches_avail = elm;
}

bool
tcache_boot(tsdn_t *tsdn)
{
	unsigned i;

	/*
	 * If necessary, clamp opt_lg_tcache_max, now that large_maxclass is
	 * known.
	 */
	if (opt_lg_tcache_max < 0 || (1U << opt_lg_tcache_max) < SMALL_MAXCLASS)
		tcache_maxclass = SMALL_MAXCLASS;
	else if ((1U << opt_lg_tcache_max) > large_maxclass)
		tcache_maxclass = large_maxclass;
	else
		tcache_maxclass = (1U << opt_lg_tcache_max);

	nhbins = size2index(tcache_maxclass) + 1;

	/* Initialize tcache_bin_info. */
	tcache_bin_info = (tcache_bin_info_t *)base_alloc(tsdn, nhbins *
	    (sizeof(tcache_bin_info_t) + sizeof(arena_cache_bin_info_t)));
	arena_cache_bin_info = (void *)tcache_bin_info +
	    sizeof(tcache_bin_info_t) * nhbins;

	if (tcache_bin_info == NULL)
		return (true);
	tcache_stack_nelms = 0;
	for (i = 0; i < NBINS; i++) {
		if ((arena_bin_info[i].nregs << 1) <= TCACHE_NSLOTS_SMALL_MIN) {
			tcache_bin_info[i].ncached_max =
			    TCACHE_NSLOTS_SMALL_MIN;
		} else if ((arena_bin_info[i].nregs << 1) <=
		    TCACHE_NSLOTS_SMALL_MAX) {
			tcache_bin_info[i].ncached_max =
			    (arena_bin_info[i].nregs << 1);
		} else {
			tcache_bin_info[i].ncached_max =
			    TCACHE_NSLOTS_SMALL_MAX;
		}
		tcache_stack_nelms += tcache_bin_info[i].ncached_max;
	}
	for (; i < nhbins; i++) {
		tcache_bin_info[i].ncached_max = TCACHE_NSLOTS_LARGE;
		tcache_stack_nelms += tcache_bin_info[i].ncached_max;
	}
	for (i = 0; i < nhbins; i++) {
		unsigned acache_max = opt_acache_size_ratio * tcache_bin_info[i].ncached_max;
		arena_cache_bin_info[i].ncached_max = acache_max;
		arena_cache_bin_info[i].flush_remain = acache_max - (acache_max >> 2);
	}

	return (false);
}
