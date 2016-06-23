/******************************************************************************/
#ifdef JEMALLOC_H_TYPES

#define	LARGE_MINCLASS		(ZU(1) << LG_LARGE_MINCLASS)

/* Maximum number of regions in one slab. */
#define	LG_SLAB_MAXREGS		(LG_PAGE - LG_TINY_MIN)
#define	SLAB_MAXREGS		(1U << LG_SLAB_MAXREGS)

/*
 * The minimum ratio of active:dirty pages per arena is computed as:
 *
 *   (nactive >> lg_dirty_mult) >= ndirty
 *
 * So, supposing that lg_dirty_mult is 3, there can be no less than 8 times as
 * many active pages as dirty pages.
 */
#define	LG_DIRTY_MULT_DEFAULT	3

#define ACACHE_EPOCH_OFF 32
#define ACACHE_EPOCH_INC ((uint64_t)1 << ACACHE_EPOCH_OFF)
#define ACACHE_LOCKBIT   ((uint64_t)1 << 31)
#define ACACHE_NCACHE_MASK (((uint64_t)1 << 31) - 1)
/* Currently fixed size. */
#define ACACHE_TCACHE_RATIO 8

typedef enum {
	purge_mode_ratio = 0,
	purge_mode_decay = 1,

	purge_mode_limit = 2
} purge_mode_t;
#define	PURGE_DEFAULT		purge_mode_ratio
/* Default decay time in seconds. */
#define	DECAY_TIME_DEFAULT	10
/* Number of event ticks between time checks. */
#define	DECAY_NTICKS_PER_UPDATE	1000

typedef struct arena_slab_data_s arena_slab_data_t;
typedef struct arena_bin_info_s arena_bin_info_t;
typedef struct arena_bin_s arena_bin_t;
typedef struct arena_s arena_t;
typedef struct arena_tdata_s arena_tdata_t;
typedef struct arena_cache_s arena_cache_t;
typedef struct arena_cache_bin_s arena_cache_bin_t;
typedef struct arena_cache_bin_info_s arena_cache_bin_info_t;


#endif /* JEMALLOC_H_TYPES */
/******************************************************************************/
#ifdef JEMALLOC_H_STRUCTS

#ifdef JEMALLOC_ARENA_STRUCTS_A
struct arena_slab_data_s {
	/* Index of bin this slab is associated with. */
	szind_t		binind;

	/* Number of free regions in slab. */
	unsigned	nfree;

	/* Per region allocated/deallocated bitmap. */
	bitmap_t	bitmap[BITMAP_GROUPS_MAX];
};
#endif /* JEMALLOC_ARENA_STRUCTS_A */

#ifdef JEMALLOC_ARENA_STRUCTS_B
/*
 * Read-only information associated with each element of arena_t's bins array
 * is stored separately, partly to reduce memory usage (only one copy, rather
 * than one per arena), but mainly to avoid false cacheline sharing.
 *
 * Each slab has the following layout:
 *
 *   /--------------------\
 *   | region 0           |
 *   |--------------------|
 *   | region 1           |
 *   |--------------------|
 *   | ...                |
 *   | ...                |
 *   | ...                |
 *   |--------------------|
 *   | region nregs-1     |
 *   \--------------------/
 */
struct arena_bin_info_s {
	/* Size of regions in a slab for this bin's size class. */
	size_t			reg_size;

	/* Total size of a slab for this bin's size class. */
	size_t			slab_size;

	/* Total number of regions in a slab for this bin's size class. */
	uint32_t		nregs;

	/*
	 * Metadata used to manipulate bitmaps for slabs associated with this
	 * bin.
	 */
	bitmap_info_t		bitmap_info;
};

struct arena_bin_s {
	/* All operations on arena_bin_t fields require lock ownership. */
	malloc_mutex_t		lock;

	/*
	 * Current slab being used to service allocations of this bin's size
	 * class.  slabcur is independent of slabs_{nonfull,full}; whenever
	 * slabcur is reassigned, the previous slab must be deallocated or
	 * inserted into slabs_{nonfull,full}.
	 */
	extent_t		*slabcur;

	/*
	 * Heap of non-full slabs.  This heap is used to assure that new
	 * allocations come from the non-full slab that is lowest in memory.
	 */
	extent_heap_t		slabs_nonfull;

	/* Ring sentinel used to track full slabs. */
	extent_t		slabs_full;

	/* Bin statistics. */
	malloc_bin_stats_t	stats;
};

struct arena_cache_bin_info_s {
	unsigned	ncached_max;	/* Upper limit on ncached. */
	unsigned	flush_remain;
};

struct arena_cache_bin_s {
	//FIXME: deal with 32bit cpus? need this in a single word.
	/*
	 * layout for data: epoch (32bits) + lock (1bit) + ncached lowest 31bits. The
	 * single bit lock is used for cache_free.
	 */
	uint64_t    data;

	/* Queued stats updated w/ atomics. */
	uint64_t    queued_nflushes;
	uint64_t    queued_nrequests;

	void		**avail;	/* Stack of available objects. */
};

struct arena_cache_s {
	/*
	 * cbins is dynamically sized (needs to be last member in this struct) and
	 * allocated contiguously off the end of arena_s.
	 */
	arena_cache_bin_t	cbins[0];
};

struct arena_s {
	/* This arena's index within the arenas array. */
	unsigned		ind;

	/*
	 * Number of threads currently assigned to this arena, synchronized via
	 * atomic operations.  Each thread has two distinct assignments, one for
	 * application-serving allocation, and the other for internal metadata
	 * allocation.  Internal metadata must not be allocated from arenas
	 * created via the arenas.extend mallctl, because the arena.<i>.reset
	 * mallctl indiscriminately discards all allocations for the affected
	 * arena.
	 *
	 *   0: Application allocation.
	 *   1: Internal metadata allocation.
	 */
	unsigned		nthreads[2];

	/*
	 * There are three classes of arena operations from a locking
	 * perspective:
	 * 1) Thread assignment (modifies nthreads) is synchronized via atomics.
	 * 2) Bin-related operations are protected by bin locks.
	 * 3) Chunk-related operations are protected by this mutex.
	 */
	malloc_mutex_t		lock;

	arena_stats_t		stats;
	/*
	 * List of tcaches for extant threads associated with this arena.
	 * Stats from these are merged incrementally, and at exit if
	 * opt_stats_print is enabled.
	 */
	ql_head(tcache_t)	tcache_ql;

	uint64_t		prof_accumbytes;

	/*
	 * PRNG state for cache index randomization of large allocation base
	 * pointers.
	 */
	uint64_t		offset_state;

	dss_prec_t		dss_prec;

	/* Minimum ratio (log base 2) of nactive:ndirty. */
	ssize_t			lg_dirty_mult;

	/* True if a thread is currently executing arena_purge_to_limit(). */
	bool			purging;

	/* Number of pages in active extents. */
	size_t			nactive;

	/*
	 * Current count of pages within unused extents that are potentially
	 * dirty, and for which madvise(... MADV_DONTNEED) has not been called.
	 * By tracking this, we can institute a limit on how much dirty unused
	 * memory is mapped for each arena.
	 */
	size_t			ndirty;

	/*
	 * Ring sentinel used to track unused dirty memory.  Dirty memory is
	 * managed as an LRU of cached extents.
	 */
	extent_t		extents_dirty;

	/*
	 * Approximate time in seconds from the creation of a set of unused
	 * dirty pages until an equivalent set of unused dirty pages is purged
	 * and/or reused.
	 */
	ssize_t			decay_time;
	/* decay_time / SMOOTHSTEP_NSTEPS. */
	nstime_t		decay_interval;
	/*
	 * Time at which the current decay interval logically started.  We do
	 * not actually advance to a new epoch until sometime after it starts
	 * because of scheduling and computation delays, and it is even possible
	 * to completely skip epochs.  In all cases, during epoch advancement we
	 * merge all relevant activity into the most recently recorded epoch.
	 */
	nstime_t		decay_epoch;
	/* decay_deadline randomness generator. */
	uint64_t		decay_jitter_state;
	/*
	 * Deadline for current epoch.  This is the sum of decay_interval and
	 * per epoch jitter which is a uniform random variable in
	 * [0..decay_interval).  Epochs always advance by precise multiples of
	 * decay_interval, but we randomize the deadline to reduce the
	 * likelihood of arenas purging in lockstep.
	 */
	nstime_t		decay_deadline;
	/*
	 * Number of dirty pages at beginning of current epoch.  During epoch
	 * advancement we use the delta between decay_ndirty and ndirty to
	 * determine how many dirty pages, if any, were generated, and record
	 * the result in decay_backlog.
	 */
	size_t			decay_ndirty;
	/*
	 * Memoized result of arena_decay_backlog_npages_limit() corresponding
	 * to the current contents of decay_backlog, i.e. the limit on how many
	 * pages are allowed to exist for the decay epochs.
	 */
	size_t			decay_backlog_npages_limit;
	/*
	 * Trailing log of how many unused dirty pages were generated during
	 * each of the past SMOOTHSTEP_NSTEPS decay epochs, where the last
	 * element is the most recent epoch.  Corresponding epoch times are
	 * relative to decay_epoch.
	 */
	size_t			decay_backlog[SMOOTHSTEP_NSTEPS];

	/* Extant large allocations. */
	ql_head(extent_t)	large;
	/* Synchronizes all large allocation/update/deallocation. */
	malloc_mutex_t		large_mtx;

	/*
	 * Heaps of extents that were previously allocated.  These are used when
	 * allocating extents, in an attempt to re-use address space.
	 */
	extent_heap_t		extents_cached[NPSIZES];
	extent_heap_t		extents_retained[NPSIZES];
	/* Protects extents_cached and extents_retained. */
	malloc_mutex_t		extents_mtx;

	/* User-configurable extent hook functions. */
	union {
		extent_hooks_t		*extent_hooks;
		void			*extent_hooks_pun;
	};

	/* Cache of extent structures that were allocated via base_alloc(). */
	ql_head(extent_t)	extent_cache;
	malloc_mutex_t		extent_cache_mtx;

	/* bins is used to store heaps of free regions. */
	arena_bin_t		bins[NBINS];
	/*
	 * Bins in acache is dynamically sized like tcache bins. acache should be
	 * the last member of this struct.
	 */
	arena_cache_t acache;
};

/* Used in conjunction with tsd for fast arena-related context lookup. */
struct arena_tdata_s {
	ticker_t		decay_ticker;
};
#endif /* JEMALLOC_ARENA_STRUCTS_B */

#endif /* JEMALLOC_H_STRUCTS */
/******************************************************************************/
#ifdef JEMALLOC_H_EXTERNS

static const size_t	large_pad =
#ifdef JEMALLOC_CACHE_OBLIVIOUS
    PAGE
#else
    0
#endif
    ;

extern purge_mode_t	opt_purge;
extern const char	*purge_mode_names[];
extern ssize_t		opt_lg_dirty_mult;
extern ssize_t		opt_decay_time;

extern const arena_bin_info_t	arena_bin_info[NBINS];
extern arena_cache_bin_info_t*arena_cache_bin_info;
extern unsigned		nhbins;

extent_t	*arena_extent_cache_alloc(tsdn_t *tsdn, arena_t *arena,
    extent_hooks_t **r_extent_hooks, void *new_addr, size_t size,
    size_t alignment, bool *zero);
void	arena_extent_cache_dalloc(tsdn_t *tsdn, arena_t *arena,
    extent_hooks_t **r_extent_hooks, extent_t *extent);
void	arena_extent_cache_maybe_insert(arena_t *arena, extent_t *extent,
    bool cache);
void	arena_extent_cache_maybe_remove(arena_t *arena, extent_t *extent,
    bool cache);
extent_t	*arena_extent_alloc_large(tsdn_t *tsdn, arena_t *arena,
    size_t usize, size_t alignment, bool *zero);
void	arena_extent_dalloc_large(tsdn_t *tsdn, arena_t *arena,
    extent_t *extent, bool locked);
void	arena_extent_ralloc_large_shrink(tsdn_t *tsdn, arena_t *arena,
    extent_t *extent, size_t oldsize);
void	arena_extent_ralloc_large_expand(tsdn_t *tsdn, arena_t *arena,
    extent_t *extent, size_t oldsize);
ssize_t	arena_lg_dirty_mult_get(tsdn_t *tsdn, arena_t *arena);
bool	arena_lg_dirty_mult_set(tsdn_t *tsdn, arena_t *arena,
    ssize_t lg_dirty_mult);
ssize_t	arena_decay_time_get(tsdn_t *tsdn, arena_t *arena);
bool	arena_decay_time_set(tsdn_t *tsdn, arena_t *arena, ssize_t decay_time);
void	arena_purge(tsdn_t *tsdn, arena_t *arena, bool all);
void	arena_maybe_purge(tsdn_t *tsdn, arena_t *arena);
void	arena_reset(tsd_t *tsd, arena_t *arena);
void	arena_tcache_fill_small(tsdn_t *tsdn, arena_t *arena,
    tcache_bin_t *tbin, szind_t binind, uint64_t prof_accumbytes);
void	arena_alloc_junk_small(void *ptr, const arena_bin_info_t *bin_info,
    bool zero);
#ifdef JEMALLOC_JET
typedef void (arena_dalloc_junk_small_t)(void *, const arena_bin_info_t *);
extern arena_dalloc_junk_small_t *arena_dalloc_junk_small;
#else
void	arena_dalloc_junk_small(void *ptr, const arena_bin_info_t *bin_info);
#endif
void	*arena_malloc_hard(tsdn_t *tsdn, arena_t *arena, size_t size,
    szind_t ind, bool zero);
void	*arena_palloc(tsdn_t *tsdn, arena_t *arena, size_t usize,
    size_t alignment, bool zero, tcache_t *tcache);
void	arena_prof_promote(tsdn_t *tsdn, extent_t *extent, const void *ptr,
    size_t usize);
void	arena_dalloc_promoted(tsdn_t *tsdn, extent_t *extent, void *ptr,
    tcache_t *tcache, bool slow_path);
void	arena_dalloc_bin_junked_locked(tsdn_t *tsdn, arena_t *arena,
    extent_t *extent, void *ptr);
void	arena_dalloc_small(tsdn_t *tsdn, arena_t *arena, extent_t *extent,
    void *ptr);
bool	arena_ralloc_no_move(tsdn_t *tsdn, extent_t *extent, void *ptr,
    size_t oldsize, size_t size, size_t extra, bool zero);
void	*arena_ralloc(tsdn_t *tsdn, arena_t *arena, extent_t *extent, void *ptr,
    size_t oldsize, size_t size, size_t alignment, bool zero, tcache_t *tcache);
dss_prec_t	arena_dss_prec_get(tsdn_t *tsdn, arena_t *arena);
bool	arena_dss_prec_set(tsdn_t *tsdn, arena_t *arena, dss_prec_t dss_prec);
ssize_t	arena_lg_dirty_mult_default_get(void);
bool	arena_lg_dirty_mult_default_set(ssize_t lg_dirty_mult);
ssize_t	arena_decay_time_default_get(void);
bool	arena_decay_time_default_set(ssize_t decay_time);
void	arena_basic_stats_merge(tsdn_t *tsdn, arena_t *arena,
    unsigned *nthreads, const char **dss, ssize_t *lg_dirty_mult,
    ssize_t *decay_time, size_t *nactive, size_t *ndirty);
void	arena_stats_merge(tsdn_t *tsdn, arena_t *arena, unsigned *nthreads,
    const char **dss, ssize_t *lg_dirty_mult, ssize_t *decay_time,
    size_t *nactive, size_t *ndirty, arena_stats_t *astats,
    malloc_bin_stats_t *bstats, malloc_large_stats_t *lstats);
unsigned	arena_nthreads_get(arena_t *arena, bool internal);
void	arena_nthreads_inc(arena_t *arena, bool internal);
void	arena_nthreads_dec(arena_t *arena, bool internal);
arena_t	*arena_new(tsdn_t *tsdn, unsigned ind);
void	arena_boot(void);
void	arena_prefork0(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork1(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork2(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork3(tsdn_t *tsdn, arena_t *arena);
void	arena_postfork_parent(tsdn_t *tsdn, arena_t *arena);
void	arena_postfork_child(tsdn_t *tsdn, arena_t *arena);

#endif /* JEMALLOC_H_EXTERNS */
/******************************************************************************/
#ifdef JEMALLOC_H_INLINES

#ifndef JEMALLOC_ENABLE_INLINE
void	arena_metadata_add(arena_t *arena, size_t size);
void	arena_metadata_sub(arena_t *arena, size_t size);
size_t	arena_metadata_get(arena_t *arena);
bool	arena_prof_accum_impl(arena_t *arena, uint64_t accumbytes);
bool	arena_prof_accum_locked(arena_t *arena, uint64_t accumbytes);
bool	arena_prof_accum(tsdn_t *tsdn, arena_t *arena, uint64_t accumbytes);
szind_t	arena_bin_index(arena_t *arena, arena_bin_t *bin);
prof_tctx_t	*arena_prof_tctx_get(tsdn_t *tsdn, const extent_t *extent,
    const void *ptr);
void	arena_prof_tctx_set(tsdn_t *tsdn, extent_t *extent, const void *ptr,
    size_t usize, prof_tctx_t *tctx);
void	arena_prof_tctx_reset(tsdn_t *tsdn, extent_t *extent, const void *ptr,
    prof_tctx_t *tctx);
void	arena_decay_ticks(tsdn_t *tsdn, arena_t *arena, unsigned nticks);
void	arena_decay_tick(tsdn_t *tsdn, arena_t *arena);
void	*arena_malloc(tsdn_t *tsdn, arena_t *arena, size_t size, szind_t ind,
    bool zero, tcache_t *tcache, bool slow_path);
arena_t	*arena_aalloc(tsdn_t *tsdn, const void *ptr);
size_t	arena_salloc(tsdn_t *tsdn, const extent_t *extent, const void *ptr);
void	arena_dalloc(tsdn_t *tsdn, extent_t *extent, void *ptr,
    tcache_t *tcache, bool slow_path);
void	arena_sdalloc(tsdn_t *tsdn, extent_t *extent, void *ptr, size_t size,
    tcache_t *tcache, bool slow_path);
void arena_cache_alloc_small(tsdn_t *tsdn, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, szind_t binind);
void *arena_cache_alloc_large(tsdn_t *tsdn, arena_t *arena, szind_t binind,
    size_t usize, bool zero);
void arena_cache_dalloc(tsdn_t *tsdn, arena_t *arena, void **items,
    size_t n_items, szind_t binind, uint64_t nrequests, const bool is_large);
void arena_cache_merge_stats(arena_t *arena, szind_t binind,
    uint64_t nrequests);
void arena_cache_flush(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
    arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
    uint64_t nrequests, const bool is_large);
#endif

#if (defined(JEMALLOC_ENABLE_INLINE) || defined(JEMALLOC_ARENA_C_))
#  ifdef JEMALLOC_ARENA_INLINE_A
JEMALLOC_INLINE void
arena_metadata_add(arena_t *arena, size_t size)
{

	atomic_add_z(&arena->stats.metadata, size);
}

JEMALLOC_INLINE void
arena_metadata_sub(arena_t *arena, size_t size)
{

	atomic_sub_z(&arena->stats.metadata, size);
}

JEMALLOC_INLINE size_t
arena_metadata_get(arena_t *arena)
{

	return (atomic_read_z(&arena->stats.metadata));
}

JEMALLOC_INLINE bool
arena_prof_accum_impl(arena_t *arena, uint64_t accumbytes)
{

	cassert(config_prof);
	assert(prof_interval != 0);

	arena->prof_accumbytes += accumbytes;
	if (arena->prof_accumbytes >= prof_interval) {
		arena->prof_accumbytes %= prof_interval;
		return (true);
	}
	return (false);
}

JEMALLOC_INLINE bool
arena_prof_accum_locked(arena_t *arena, uint64_t accumbytes)
{

	cassert(config_prof);

	if (likely(prof_interval == 0))
		return (false);
	return (arena_prof_accum_impl(arena, accumbytes));
}

JEMALLOC_INLINE bool
arena_prof_accum(tsdn_t *tsdn, arena_t *arena, uint64_t accumbytes)
{

	cassert(config_prof);

	if (likely(prof_interval == 0))
		return (false);

	{
		bool ret;

		malloc_mutex_lock(tsdn, &arena->lock);
		ret = arena_prof_accum_impl(arena, accumbytes);
		malloc_mutex_unlock(tsdn, &arena->lock);
		return (ret);
	}
}
#  endif /* JEMALLOC_ARENA_INLINE_A */

#  ifdef JEMALLOC_ARENA_INLINE_B
JEMALLOC_INLINE szind_t
arena_bin_index(arena_t *arena, arena_bin_t *bin)
{
	szind_t binind = (szind_t)(bin - arena->bins);
	assert(binind < NBINS);
	return (binind);
}

JEMALLOC_INLINE prof_tctx_t *
arena_prof_tctx_get(tsdn_t *tsdn, const extent_t *extent, const void *ptr)
{

	cassert(config_prof);
	assert(ptr != NULL);

	if (unlikely(!extent_slab_get(extent)))
		return (large_prof_tctx_get(tsdn, extent));
	return ((prof_tctx_t *)(uintptr_t)1U);
}

JEMALLOC_INLINE void
arena_prof_tctx_set(tsdn_t *tsdn, extent_t *extent, const void *ptr,
    size_t usize, prof_tctx_t *tctx)
{

	cassert(config_prof);
	assert(ptr != NULL);

	if (unlikely(!extent_slab_get(extent)))
		large_prof_tctx_set(tsdn, extent, tctx);
}

JEMALLOC_INLINE void
arena_prof_tctx_reset(tsdn_t *tsdn, extent_t *extent, const void *ptr,
    prof_tctx_t *tctx)
{

	cassert(config_prof);
	assert(ptr != NULL);
	assert(!extent_slab_get(extent));

	large_prof_tctx_reset(tsdn, extent);
}

JEMALLOC_ALWAYS_INLINE void
arena_decay_ticks(tsdn_t *tsdn, arena_t *arena, unsigned nticks)
{
	tsd_t *tsd;
	ticker_t *decay_ticker;

	if (unlikely(tsdn_null(tsdn)))
		return;
	tsd = tsdn_tsd(tsdn);
	decay_ticker = decay_ticker_get(tsd, arena->ind);
	if (unlikely(decay_ticker == NULL))
		return;
	if (unlikely(ticker_ticks(decay_ticker, nticks)))
		arena_purge(tsdn, arena, false);
}

JEMALLOC_ALWAYS_INLINE void
arena_decay_tick(tsdn_t *tsdn, arena_t *arena)
{

	arena_decay_ticks(tsdn, arena, 1);
}

JEMALLOC_ALWAYS_INLINE void *
arena_malloc(tsdn_t *tsdn, arena_t *arena, size_t size, szind_t ind, bool zero,
    tcache_t *tcache, bool slow_path)
{

	assert(!tsdn_null(tsdn) || tcache == NULL);
	assert(size != 0);

	if (likely(tcache != NULL)) {
		if (likely(size <= SMALL_MAXCLASS)) {
			return (tcache_alloc_small(tsdn_tsd(tsdn), arena,
			    tcache, size, ind, zero, slow_path));
		}
		if (likely(size <= tcache_maxclass)) {
			return (tcache_alloc_large(tsdn_tsd(tsdn), arena,
			    tcache, size, ind, zero, slow_path));
		}
		/* (size > tcache_maxclass) case falls through. */
		assert(size > tcache_maxclass);
	}

	return (arena_malloc_hard(tsdn, arena, size, ind, zero));
}

JEMALLOC_ALWAYS_INLINE arena_t *
arena_aalloc(tsdn_t *tsdn, const void *ptr)
{

	return (extent_arena_get(iealloc(tsdn, ptr)));
}

/* Return the size of the allocation pointed to by ptr. */
JEMALLOC_ALWAYS_INLINE size_t
arena_salloc(tsdn_t *tsdn, const extent_t *extent, const void *ptr)
{
	size_t ret;

	assert(ptr != NULL);

	if (likely(extent_slab_get(extent)))
		ret = index2size(extent_slab_data_get_const(extent)->binind);
	else
		ret = large_salloc(tsdn, extent);

	return (ret);
}

JEMALLOC_ALWAYS_INLINE void
arena_dalloc(tsdn_t *tsdn, extent_t *extent, void *ptr, tcache_t *tcache,
    bool slow_path)
{

	assert(!tsdn_null(tsdn) || tcache == NULL);
	assert(ptr != NULL);

	if (likely(extent_slab_get(extent))) {
		/* Small allocation. */
		if (likely(tcache != NULL)) {
			szind_t binind = extent_slab_data_get(extent)->binind;
			tcache_dalloc_small(tsdn_tsd(tsdn), tcache, ptr, binind,
			    slow_path);
		} else {
			arena_dalloc_small(tsdn, extent_arena_get(extent),
			    extent, ptr);
		}
	} else {
		size_t usize = extent_usize_get(extent);

		if (likely(tcache != NULL) && usize <= tcache_maxclass) {
			if (config_prof && unlikely(usize <= SMALL_MAXCLASS)) {
				arena_dalloc_promoted(tsdn, extent, ptr,
				    tcache, slow_path);
			} else {
				tcache_dalloc_large(tsdn_tsd(tsdn), tcache,
				    ptr, usize, slow_path);
			}
		} else
			large_dalloc(tsdn, extent);
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_sdalloc(tsdn_t *tsdn, extent_t *extent, void *ptr, size_t size,
    tcache_t *tcache, bool slow_path)
{

	assert(!tsdn_null(tsdn) || tcache == NULL);
	assert(ptr != NULL);

	if (likely(extent_slab_get(extent))) {
		/* Small allocation. */
		if (likely(tcache != NULL)) {
			szind_t binind = size2index(size);
			assert(binind == extent_slab_data_get(extent)->binind);
			tcache_dalloc_small(tsdn_tsd(tsdn), tcache, ptr, binind,
			    slow_path);
		} else {
			arena_dalloc_small(tsdn, extent_arena_get(extent),
			    extent, ptr);
		}
	} else {
		if (likely(tcache != NULL) && size <= tcache_maxclass) {
			if (config_prof && unlikely(size <= SMALL_MAXCLASS)) {
				arena_dalloc_promoted(tsdn, extent, ptr,
				    tcache, slow_path);
			} else {
				tcache_dalloc_large(tsdn_tsd(tsdn), tcache, ptr,
				    size, slow_path);
			}
		} else
			large_dalloc(tsdn, extent);
	}
}

JEMALLOC_ALWAYS_INLINE uint64_t
cbin_info_get(arena_cache_bin_t *cbin, size_t *ncached, bool *locked)
{
	uint64_t val = ACCESS_ONCE(cbin->data);

	*ncached = val & ACACHE_NCACHE_MASK;
	if (val & ACACHE_LOCKBIT)
		*locked = true;
	else
		*locked = false;

	return val;
}

/* Return number of items allocated */
JEMALLOC_ALWAYS_INLINE size_t
cbin_alloc_to_tbin(arena_cache_bin_t *cbin, tcache_bin_t *tbin, unsigned n_alloc)
{
	uint64_t val, new_val;
	size_t ncached;
	bool cas_fail, locked;
retry:
	val = cbin_info_get(cbin, &ncached, &locked);
	if (locked || (ncached == 0))
		return 0;

	if (n_alloc > ncached)
		n_alloc = ncached;
	/* tbin->avail[-ncached ... -1] are available items. */
	memcpy(&tbin->avail[-(int)n_alloc], &cbin->avail[ncached - n_alloc],
	    sizeof(void *) * n_alloc);

	/* Updated value: epoch and ncached */
	new_val = val + ACACHE_EPOCH_INC - n_alloc;
	cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
	if (unlikely(cas_fail)) {
		assert((ACCESS_ONCE(cbin->data) >> ACACHE_EPOCH_OFF) !=
		    (val >> ACACHE_EPOCH_OFF));
		goto retry;
	}
	assert(tbin->avail[-(int)n_alloc] && tbin->avail[-1]);

	return n_alloc;
}

/*
 * Allocate a single item from acache if available. Currently used for large
 * allocations only.
 */
JEMALLOC_ALWAYS_INLINE void *
cbin_alloc(arena_cache_bin_t *cbin)
{
	uint64_t val, new_val;
	size_t ncached;
	bool cas_fail, locked;
	void *ret;
retry:
	val = cbin_info_get(cbin, &ncached, &locked);
	if (locked || (ncached == 0))
		return 0;

	ret = cbin->avail[ncached - 1];
	/* Updated value: epoch and ncache */
	new_val = val + ACACHE_EPOCH_INC - 1;
	cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
	if (unlikely(cas_fail)) {
		assert((ACCESS_ONCE(cbin->data) >> ACACHE_EPOCH_OFF) !=
		    (val >> ACACHE_EPOCH_OFF));
		goto retry;
	}

	return ret;
}

/* Allocated a single item. Use cached memory if available. */
JEMALLOC_ALWAYS_INLINE void *
arena_cache_alloc_large(tsdn_t *tsdn, arena_t *arena, szind_t binind,
    size_t usize, bool zero)
{
	arena_cache_bin_t *cbin = &arena->acache.cbins[binind];
	void *ret = NULL;

	assert(binind <= nhbins);
	if (!zero)
		ret = cbin_alloc(cbin);

	return ret ? ret : large_malloc(tsdn, arena, usize, zero);
}

/*
 * Fill tcache from acache if there is memory cached in acache. Otherwise, fill
 * from arena directly.
 */
JEMALLOC_ALWAYS_INLINE void
arena_cache_alloc_small(tsdn_t *tsdn, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, szind_t binind)
{
	arena_cache_bin_t *cbin = &arena->acache.cbins[binind];
	size_t allocated, nfill;

	assert(binind < NBINS);
	nfill = tcache_bin_info[binind].ncached_max >> tbin->lg_fill_div;
	assert(nfill);
	allocated = cbin_alloc_to_tbin(cbin, tbin, nfill);
	if (allocated) {
		tbin->ncached = allocated;
		return;
	}

	/* fill through arena directly */
	arena_tcache_fill_small(tsdn, arena, tbin, binind, config_prof ?
	    tcache->prof_accumbytes : 0);
	if (config_prof)
		tcache->prof_accumbytes = 0;
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_merge_stats(arena_t *arena, szind_t binind, uint64_t nrequests)
{
	arena_cache_bin_t *cbin = &arena->acache.cbins[binind];
	/* To reduce mutex contention, delay the stats updates. */
	atomic_add_uint64(&cbin->queued_nrequests, nrequests);
	if (binind < NBINS)
		atomic_add_uint64(&cbin->queued_nflushes, 1);
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_stats_merge_locked(arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, szind_t binind, uint64_t nrequests,
		const bool is_large)
{
	uint64_t queued_nrequests = ACCESS_ONCE(cbin->queued_nrequests);

	if (nrequests || queued_nrequests) {
		atomic_sub_uint64(&cbin->queued_nrequests, queued_nrequests);

		if (is_large) {
			arena->stats.nrequests_large += nrequests + queued_nrequests;
			arena->stats.lstats[binind - NBINS].nrequests +=
				nrequests + queued_nrequests;
		} else {
			uint64_t queued_nflushes = ACCESS_ONCE(cbin->queued_nflushes);
			atomic_sub_uint64(&cbin->queued_nflushes, queued_nflushes);
			bin->stats.nflushes += 1 + queued_nflushes;
			bin->stats.nrequests += nrequests + queued_nrequests;
		}
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_flush(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
		uint64_t nrequests, const bool is_large)
{

	if (is_large)
		malloc_mutex_lock(tsdn, &arena->lock);
	else
		malloc_mutex_lock(tsdn, &bin->lock);
	if (config_stats) {
		arena_cache_stats_merge_locked(arena, bin, cbin, binind, nrequests,
		    is_large);
	}

	for (size_t i = 0; i < n_items; i++) {
		void *ptr = items[i];
		extent_t *extent = iealloc(tsdn, ptr);

		if (is_large) {
			large_dalloc_junked_locked(tsdn, extent);
		} else {
			arena_dalloc_bin_junked_locked(tsdn, arena, extent, ptr);
		}
	}
	if (is_large)
		malloc_mutex_unlock(tsdn, &arena->lock);
	else
		malloc_mutex_unlock(tsdn, &bin->lock);
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_dalloc(tsdn_t *tsdn, arena_t *arena, void **items,
    size_t n_items, szind_t binind, uint64_t nrequests, const bool is_large)
{
	arena_bin_t *bin = &arena->bins[binind];
	arena_cache_bin_t *cbin = &arena->acache.cbins[binind];
	uint64_t val, new_val;
	size_t ncached;
	bool cas_fail, bin_locked = false;

retry:
	val = cbin_info_get(cbin, &ncached, &bin_locked);
	if (bin_locked) {
		/* Flush back to arena if acache is locked. */
		arena_cache_flush(tsdn, arena, bin, cbin, items, n_items, binind,
		    nrequests, is_large);
		return;
	}

	new_val = (val + ACACHE_EPOCH_INC) | ACACHE_LOCKBIT;
	assert(n_items <= arena_cache_bin_info[binind].ncached_max);
	if (ncached + n_items > arena_cache_bin_info[binind].ncached_max) {
		/* We need to flush some items back to arena. */
		size_t nkeep;
		if (n_items >= arena_cache_bin_info[binind].flush_remain) {
			nkeep = 0;
		} else {
			nkeep = arena_cache_bin_info[binind].flush_remain - n_items;
		}
		/* After flushing back to arena, and return the n_items to acache, we leave
		 * flush_remain items in acache. */
		assert(nkeep < ncached);

		/*
		 * If we use cbin_alloc to get items from the cbin, it's possible to avoid
		 * locking the cbin here. However it requires additional memcpy.
		 */
		/* lock the cache bin */
		cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
		if (unlikely(cas_fail)) {
			goto retry;
		}
		bin_locked = true;
		arena_cache_flush(tsdn, arena, bin, cbin, &cbin->avail[nkeep],
		    ncached - nkeep, binind, nrequests, is_large);

		/* Update new_val which will be used for unlocking. */
		new_val -= ncached - nkeep;
		ncached = nkeep;
		/* Fall through to free to cbin and unlock. */
	}

	assert(ncached + n_items <= arena_cache_bin_info[binind].ncached_max);
	if (!bin_locked) {
		/* lock the cache bin */
		cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
		if (unlikely(cas_fail)) {
			goto retry;
		}
	}

	/*
	 * Do the memcpy, then release the bitlock & update ncached in a single
	 * store.
	 */
	memcpy(&cbin->avail[ncached], items, sizeof(void *) * n_items);
	cc_barrier(); /* For x86, compiler barrier is sufficient. */
	/* Release the bitlock, and increment ncached. */
	assert(cbin->avail[ncached + n_items - 1]);
	cbin->data = (new_val & ~ACACHE_LOCKBIT) + n_items;
}

#  endif /* JEMALLOC_ARENA_INLINE_B */
#endif

#endif /* JEMALLOC_H_INLINES */
/******************************************************************************/
