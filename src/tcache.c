#define	JEMALLOC_TCACHE_C_
#include "jemalloc/internal/jemalloc_internal.h"

/******************************************************************************/
/* Data. */

bool	opt_tcache = true;
ssize_t	opt_lg_tcache_max = LG_TCACHE_MAXCLASS_DEFAULT;

tcache_bin_info_t	*tcache_bin_info;
static unsigned		stack_nelms; /* Total stack elms per tcache. */

unsigned		nhbins;
size_t			tcache_maxclass;

tcaches_t		*tcaches;
ccache_t		**ccaches;

/* Index of first element within tcaches that has never been used. */
static unsigned		tcaches_past;

/* Head of singly linked list tracking available tcaches elements. */
static tcaches_t	*tcaches_avail;

/******************************************************************************/

size_t	tcache_salloc(const void *ptr)
{

	return (arena_salloc(ptr, false));
}

volatile uint64_t nevents;

void
tcache_event_hard(tsd_t *tsd, tcache_t *tcache)
{
	szind_t binind;
    ccache_bin_t *cbin;
	tcache_bin_t *tbin;
	tcache_bin_info_t *tbin_info;
    uint32_t low_water, ncached, tid;
    uint64_t val;

    ccache_t *ccache = tsd_ccache_get(tsd);

    /* TODO: fix this ugly hack. */
    binind = ccache->next_gc_bin;
    if (binind >= nhbins) {
        ccache->next_gc_bin = 0;
        return;
    }

    cbin = &(tsd_ccache_get(tsd)->cbins[binind]);
	tbin = &tcache->tbins[binind];
	tbin_info = &tcache_bin_info[binind];

    ccache_bin_lock(tsd, cbin);

    val = ACCESS_ONCE(cbin->data);
    cache_bin_get_info(val, &ncached, &low_water, &tid);

    /* if (binind == 30) { */
    /*     __sync_fetch_and_add(&nevents, 1); */
    /*     printf(">>>>>>>>>>> thd %p, now %llu events: low_water %d (refilled %d), ncached %d div %d\n", */
    /*            tsd, nevents, low_water, cbin->refilled, ncached, cbin->lg_fill_div); */
    /* } */

	if (low_water > 0) {
		/*
		 * Flush (ceiling) 3/4 of the objects below the low water mark.
		 */
		if (binind < NBINS) {
			val = tcache_bin_flush_small_locked(tsd, tcache, tbin, cbin, binind,
			    ncached - low_water + (low_water >> 2));
		} else {
			val = tcache_bin_flush_large_locked(tsd, tcache, tbin, cbin, binind,
                ncached - low_water + (low_water >> 2));
		}
        /* Get the updated ncached value. */
        cache_bin_get_info(val, &ncached, &low_water, &tid);

		/*
		 * Reduce fill count by 2X.  Limit lg_fill_div such that the
		 * fill count is always at least 1.
		 */
		if ((tbin_info->ncached_max >> (cbin->lg_fill_div+1)) >= 1) {
			cbin->lg_fill_div++;
//            if (binind == 30) printf("div ++ : %d\n", cbin->lg_fill_div);
        }
	} else if (low_water == 0 && cbin->refilled) {
		/*
		 * Increase fill count by 2X.  Make sure lg_fill_div stays
		 * greater than 0.
		 */
		if (cbin->lg_fill_div > 1) {
			cbin->lg_fill_div--;
//            if (binind == 30) printf("div ------- : %d\n", cbin->lg_fill_div);
        }
	} else {
//        if (binind == 30) printf("nothing: div %d, low %d, refilled %d\n", cbin->lg_fill_div, low_water, cbin->refilled);
    }

    cbin->refilled = false;
    /* Update low_water to ncached. */
    val = cbin_data_pack(ncached, ncached, tid);
    assert(!CCACHE_IS_LOCKED(val));
    ccache_bin_unlock(tsd, cbin, val);

	ccache->next_gc_bin++;
	if (ccache->next_gc_bin >= nhbins)
		ccache->next_gc_bin = 0;

    percpu_cache_arena_update(tsd, malloc_getcpu());
}

void *
tcache_alloc_small_hard(tsd_t *tsd, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, ccache_bin_t *cbin, szind_t binind, uint32_t curr_thd,
    ret_status_t *ret_status)
{
	void *ret;

	arena_tcache_fill_small(tsd, arena, tbin, cbin, binind, config_prof ?
	    tcache->prof_accumbytes : 0);
	if (config_prof)
		tcache->prof_accumbytes = 0;
	ret = tcache_alloc_easy(cbin, curr_thd, ret_status);

	return (ret);
}

uint64_t
tcache_bin_flush_small_locked(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem)
{
	arena_t *arena;
	void *ptr;
	unsigned i, nflush, ndeferred;
	bool merged_stats = false;
    uint32_t low_water, ncached, tid;
    uint64_t val;

    val = ACCESS_ONCE(cbin->data);
    assert(val & CBIN_BITLOCK);
    cache_bin_get_info(val, &ncached, &low_water, &tid);
//ifdef percpu
    if (rem >= ncached) {
        val &= ~CBIN_BITLOCK;
        goto label_return;
    }
//elseif
//	assert(rem <= tbin->ncached);

	arena = arena_choose(tsd, NULL);
	assert(arena != NULL);
	for (nflush = ncached - rem; nflush > 0; nflush = ndeferred) {
		/* Lock the arena bin associated with the first object. */
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(
		    *(cbin->avail - 1));
		arena_t *bin_arena = extent_node_arena_get(&chunk->node);
		arena_bin_t *bin = &bin_arena->bins[binind];

		if (config_prof && bin_arena == arena) {
			if (arena_prof_accum(arena, tcache->prof_accumbytes))
				prof_idump();
			tcache->prof_accumbytes = 0;
		}

		malloc_mutex_lock(&bin->lock);
		if (config_stats && bin_arena == arena) {
			assert(!merged_stats);
			merged_stats = true;
			bin->stats.nflushes++;
			bin->stats.nrequests += tbin->tstats.nrequests;
			tbin->tstats.nrequests = 0;
		}
		ndeferred = 0;
		for (i = 0; i < nflush; i++) {
			ptr = *(cbin->avail - 1 - i);
			assert(ptr != NULL);
			chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
			if (extent_node_arena_get(&chunk->node) == bin_arena) {
				size_t pageind = ((uintptr_t)ptr -
				    (uintptr_t)chunk) >> LG_PAGE;
				arena_chunk_map_bits_t *bitselm =
				    arena_bitselm_get(chunk, pageind);
				arena_dalloc_bin_junked_locked(bin_arena, chunk,
				    ptr, bitselm);
			} else {
				/*
				 * This object was allocated via a different
				 * arena bin than the one that is currently
				 * locked.  Stash the object, so that it can be
				 * handled in a future pass.
				 */
				*(cbin->avail - 1 - ndeferred) = ptr;
				ndeferred++;
			}
		}
		malloc_mutex_unlock(&bin->lock);
		arena_decay_ticks(tsd, bin_arena, nflush - ndeferred);
	}
	if (config_stats && !merged_stats) {
		/*
		 * The flush loop didn't happen to flush to this thread's
		 * arena, so the stats didn't get merged.  Manually do so now.
		 */
		arena_bin_t *bin = &arena->bins[binind];
		malloc_mutex_lock(&bin->lock);
		bin->stats.nflushes++;
		bin->stats.nrequests += tbin->tstats.nrequests;
		tbin->tstats.nrequests = 0;
		malloc_mutex_unlock(&bin->lock);
	}

	memmove(cbin->avail - rem, cbin->avail - ncached, rem *
	    sizeof(void *));

    ncached = rem;
    if (ncached < low_water)
		low_water = ncached;
    val = cbin_data_pack(rem, low_water, tid);
label_return:
    return val;
}

void
tcache_bin_flush_small(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem)
{
    uint64_t bin_data;

    ccache_bin_lock(tsd, cbin);
    bin_data = tcache_bin_flush_small_locked(tsd, tcache, tbin, cbin, binind, rem);
    ccache_bin_unlock(tsd, cbin, bin_data);
}

uint64_t
tcache_bin_flush_large_locked(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem)
{
	arena_t *arena;
	void *ptr;
	unsigned i, nflush, ndeferred;
    uint32_t low_water, ncached, tid;
    uint64_t val;
	bool merged_stats = false;

    val = ACCESS_ONCE(cbin->data);
    cache_bin_get_info(val, &ncached, &low_water, &tid);
//ifdef percpu
    if (rem >= ncached) {
        val &= ~CBIN_BITLOCK;
        goto label_return;
    }
//elseif
//	assert(rem <= tbin->ncached);

	assert(binind < nhbins);
	arena = arena_choose(tsd, NULL);
	assert(arena != NULL);
	for (nflush = ncached - rem; nflush > 0; nflush = ndeferred) {
		/* Lock the arena associated with the first object. */
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(
		    *(cbin->avail - 1));
		arena_t *locked_arena = extent_node_arena_get(&chunk->node);
		UNUSED bool idump;

		if (config_prof)
			idump = false;
		malloc_mutex_lock(&locked_arena->lock);
		if ((config_prof || config_stats) && locked_arena == arena) {
			if (config_prof) {
				idump = arena_prof_accum_locked(arena,
				    tcache->prof_accumbytes);
				tcache->prof_accumbytes = 0;
			}
			if (config_stats) {
				merged_stats = true;
				arena->stats.nrequests_large +=
				    tbin->tstats.nrequests;
				arena->stats.lstats[binind - NBINS].nrequests +=
				    tbin->tstats.nrequests;
				tbin->tstats.nrequests = 0;
			}
		}
		ndeferred = 0;
		for (i = 0; i < nflush; i++) {
			ptr = *(cbin->avail - 1 - i);
			assert(ptr != NULL);
			chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
			if (extent_node_arena_get(&chunk->node) ==
			    locked_arena) {
				arena_dalloc_large_junked_locked(locked_arena,
				    chunk, ptr);
			} else {
				/*
				 * This object was allocated via a different
				 * arena than the one that is currently locked.
				 * Stash the object, so that it can be handled
				 * in a future pass.
				 */
				*(cbin->avail - 1 - ndeferred) = ptr;
				ndeferred++;
			}
		}
		malloc_mutex_unlock(&locked_arena->lock);
		if (config_prof && idump)
			prof_idump();
		arena_decay_ticks(tsd, locked_arena, nflush - ndeferred);
	}
	if (config_stats && !merged_stats) {
		/*
		 * The flush loop didn't happen to flush to this thread's
		 * arena, so the stats didn't get merged.  Manually do so now.
		 */
		malloc_mutex_lock(&arena->lock);
		arena->stats.nrequests_large += tbin->tstats.nrequests;
		arena->stats.lstats[binind - NBINS].nrequests +=
		    tbin->tstats.nrequests;
		tbin->tstats.nrequests = 0;
		malloc_mutex_unlock(&arena->lock);
	}

	memmove(cbin->avail - rem, cbin->avail - ncached, rem *
	    sizeof(void *));

    ncached = rem;
    if (ncached < low_water)
		low_water = ncached;
    val = cbin_data_pack(rem, low_water, tid);
label_return:
    return val;
}

void
tcache_bin_flush_large(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem)
{
    uint64_t bin_data;

    ccache_bin_lock(tsd, cbin);
    bin_data = tcache_bin_flush_large_locked(tsd, tcache, tbin, cbin, binind, rem);
    ccache_bin_unlock(tsd, cbin, bin_data);
}

void
tcache_arena_associate(tcache_t *tcache, arena_t *arena)
{

	if (config_stats) {
		/* Link into list of extant tcaches. */
		malloc_mutex_lock(&arena->lock);
		ql_elm_new(tcache, link);
		ql_tail_insert(&arena->tcache_ql, tcache, link);
		malloc_mutex_unlock(&arena->lock);
	}
}

void
tcache_arena_reassociate(tcache_t *tcache, arena_t *oldarena, arena_t *newarena)
{

	tcache_arena_dissociate(tcache, oldarena);
	tcache_arena_associate(tcache, newarena);
}

void
tcache_arena_dissociate(tcache_t *tcache, arena_t *arena)
{

	if (config_stats) {
		/* Unlink from list of extant tcaches. */
		malloc_mutex_lock(&arena->lock);
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
		tcache_stats_merge(tcache, arena);
		malloc_mutex_unlock(&arena->lock);
	}
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
	return (tcache_create(tsd, arena));
}

#include <sys/syscall.h>
tcache_t *
tcache_create(tsd_t *tsd, arena_t *arena)
{
	tcache_t *tcache;
	size_t size/* , stack_offset */;

	size = offsetof(tcache_t, tbins) + (sizeof(tcache_bin_t) * nhbins);
	/* Naturally align the pointer stacks. */
	size = PTR_CEILING(size);

    // TODO only for tcache
	/* stack_offset = size; */
	/* size += stack_nelms * sizeof(void *); */

	/* Avoid false cacheline sharing (additional cacheline for prefetching). */
	size = sa2u(size + CACHELINE, CACHELINE);

	tcache = ipallocztm(tsd, size, CACHELINE, true, false, true,
	    arena_get(0, false));
	if (tcache == NULL)
		return (NULL);

	tcache_arena_associate(tcache, arena);

	ticker_init(&tcache->gc_ticker, TCACHE_GC_INCR);

	assert((TCACHE_NSLOTS_SMALL_MAX & 1U) == 0);
//#ifndef
	/* for (i = 0; i < nhbins; i++) { */
	/* 	stack_offset += tcache_bin_info[i].ncached_max * sizeof(void *); */
	/* 	/\* */
	/* 	 * avail points past the available space.  Allocations will */
	/* 	 * access the slots toward higher addresses (for the benefit of */
	/* 	 * prefetch). */
	/* 	 *\/ */
	/* 	tcache->tbins[i].avail = (void **)((uintptr_t)tcache + */
	/* 	    (uintptr_t)stack_offset); */
	/* } */

	return (tcache);
}

volatile unsigned t_killed = 0;

static void
tcache_destroy(tsd_t *tsd, tcache_t *tcache)
{
	arena_t *arena;
//    ccache_t *ccache;
	unsigned i;

	arena = arena_choose(tsd, NULL);
	tcache_arena_dissociate(tcache, arena);

    __sync_fetch_and_add(&t_killed, 1);

//    ccache = tsd_ccache_get(tsd);

	for (i = 0; i < NBINS; i++) {
		tcache_bin_t *tbin = &tcache->tbins[i];
//        ccache_bin_t *cbin = &ccache->cbins[i];

//        tcache_bin_flush_small(tsd, tcache, tbin, cbin, i, 0);
		if (config_stats && tbin->tstats.nrequests != 0) {
			arena_bin_t *bin = &arena->bins[i];
			malloc_mutex_lock(&bin->lock);
			bin->stats.nrequests += tbin->tstats.nrequests;
			malloc_mutex_unlock(&bin->lock);
		}
	}

	for (; i < nhbins; i++) {
		tcache_bin_t *tbin = &tcache->tbins[i];
//        ccache_bin_t *cbin = &ccache->cbins[i];

//        tcache_bin_flush_large(tsd, tcache, tbin, cbin, i, 0);
		if (config_stats && tbin->tstats.nrequests != 0) {
			malloc_mutex_lock(&arena->lock);
			arena->stats.nrequests_large += tbin->tstats.nrequests;
			arena->stats.lstats[i - NBINS].nrequests +=
			    tbin->tstats.nrequests;
			malloc_mutex_unlock(&arena->lock);
		}
	}

	if (config_prof && tcache->prof_accumbytes > 0 &&
	    arena_prof_accum(arena, tcache->prof_accumbytes))
		prof_idump();

	idalloctm(tsd, tcache, false, true, true);
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

/* Caller must own arena->lock. */
void
tcache_stats_merge(tcache_t *tcache, arena_t *arena)
{
	unsigned i;

	cassert(config_stats);

	/* Merge and reset tcache stats. */
	for (i = 0; i < NBINS; i++) {
		arena_bin_t *bin = &arena->bins[i];
		tcache_bin_t *tbin = &tcache->tbins[i];
		malloc_mutex_lock(&bin->lock);
		bin->stats.nrequests += tbin->tstats.nrequests;
		malloc_mutex_unlock(&bin->lock);
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
tcaches_create(tsd_t *tsd, unsigned *r_ind)
{
	tcache_t *tcache;
	tcaches_t *elm;

	if (tcaches == NULL) {
		tcaches = base_alloc(sizeof(tcache_t *) *
		    (MALLOCX_TCACHE_MAX+1));
		if (tcaches == NULL)
			return (true);
	}

	if (tcaches_avail == NULL && tcaches_past > MALLOCX_TCACHE_MAX)
		return (true);
	tcache = tcache_create(tsd, arena_get(0, false));
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

ccache_t *
ccache_create(tsd_t *tsd)
{
	ccache_t *ccache;
	size_t size, stack_offset;
	unsigned i;

	size = offsetof(ccache_t, cbins) + sizeof(ccache_bin_t) * nhbins;
	/* Naturally align the pointer stacks. */
	size = PTR_CEILING(size);
	stack_offset = size;
	size += stack_nelms * sizeof(void *);
	/* Avoid false cacheline sharing. 1 extra cacheline for prefeteching. */
	size = sa2u(size + CACHELINE, CACHELINE);

	ccache = ipallocztm(tsd, size, CACHELINE, true, false, true,
	    arena_get(0, false));
	if (ccache == NULL)
		return (NULL);

    ccache->next_gc_bin = 0;
	assert((TCACHE_NSLOTS_SMALL_MAX & 1U) == 0);
	for (i = 0; i < nhbins; i++) {
		ccache->cbins[i].lg_fill_div = 1;
		ccache->cbins[i].refilled = false;
		stack_offset += tcache_bin_info[i].ncached_max * sizeof(void *);
		/*
		 * avail points past the available space.  Allocations will
		 * access the slots toward higher addresses (for the benefit of
		 * prefetch).
		 */
		ccache->cbins[i].avail = (void **)((uintptr_t)ccache +
		    (uintptr_t)stack_offset);
	}

	return (ccache);
}

bool
ccaches_create(unsigned ncaches)
{
	ccache_t *ccache;
    unsigned i;
    tsd_t *tsd = tsd_fetch();

    /* ccaches is only accessed when thread migration is detected. */
    ccaches = base_alloc(sizeof(ccache_t *) * ncaches);
    if (ccaches == NULL)
        return (true);

    for (i = 0; i < ncaches; i++) {
        ccache = ccache_create(tsd);
        if (ccache == NULL)
            return (true);
        ccache->ind = i;
        ccaches[i] = ccache;
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
tcache_boot(void)
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
	tcache_bin_info = (tcache_bin_info_t *)base_alloc(nhbins *
	    sizeof(tcache_bin_info_t));
	if (tcache_bin_info == NULL)
		return (true);
	stack_nelms = 0;
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
		stack_nelms += tcache_bin_info[i].ncached_max;
	}
	for (; i < nhbins; i++) {
		tcache_bin_info[i].ncached_max = TCACHE_NSLOTS_LARGE;
		stack_nelms += tcache_bin_info[i].ncached_max;
	}

	return (false);
}
