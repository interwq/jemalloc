/******************************************************************************/
#ifdef JEMALLOC_H_TYPES

#define	LARGE_MINCLASS		(ZU(1) << LG_LARGE_MINCLASS)

/* Maximum number of regions in one run. */
#define	LG_RUN_MAXREGS		(LG_PAGE - LG_TINY_MIN)
#define	RUN_MAXREGS		(1U << LG_RUN_MAXREGS)

/*
 * Minimum redzone size.  Redzones may be larger than this if necessary to
 * preserve region alignment.
 */
#define	REDZONE_MINSIZE		16

/*
 * The minimum ratio of active:dirty pages per arena is computed as:
 *
 *   (nactive >> lg_dirty_mult) >= ndirty
 *
 * So, supposing that lg_dirty_mult is 3, there can be no less than 8 times as
 * many active pages as dirty pages.
 */
#define	LG_DIRTY_MULT_DEFAULT	3

typedef enum {
	purge_mode_ratio = 0,
	purge_mode_decay = 1,

	purge_mode_limit = 2
} purge_mode_t;
#define	PURGE_DEFAULT		purge_mode_decay
/* Default decay time in seconds. */
#define	DECAY_TIME_DEFAULT	10
/* Number of event ticks between time checks. */
#define	DECAY_NTICKS_PER_UPDATE	1000

typedef enum {
	percpu_arena_disable = 0,
	percpu_arena_enable = 1,
	per_phycpu_arena_enable = 2, 	/* i.e. hyper threads share arena */

	percpu_arena_mode_limit = 3
} percpu_arena_mode_t;
#define PERCPU_ARENA_DEFAULT percpu_arena_disable
#define PURGE_THREAD_DEFAULT false
/* # of seconds between purging. */
#define PURGE_THREAD_INTERVAL 1

typedef struct arena_runs_dirty_link_s arena_runs_dirty_link_t;
typedef struct arena_avail_links_s arena_avail_links_t;
typedef struct arena_run_s arena_run_t;
typedef struct arena_chunk_map_bits_s arena_chunk_map_bits_t;
typedef struct arena_chunk_map_misc_s arena_chunk_map_misc_t;
typedef struct arena_chunk_s arena_chunk_t;
typedef struct arena_bin_info_s arena_bin_info_t;
typedef struct arena_bin_s arena_bin_t;
typedef struct arena_s arena_t;
typedef struct arena_tdata_s arena_tdata_t;
typedef struct arena_cache_s arena_cache_t;
typedef struct arena_cache_bin_s arena_cache_bin_t;
typedef struct arena_cache_bin_info_s arena_cache_bin_info_t;
typedef uint64_t acache_state_t;


#endif /* JEMALLOC_H_TYPES */
/******************************************************************************/
#ifdef JEMALLOC_H_STRUCTS

#ifdef JEMALLOC_ARENA_STRUCTS_A
struct arena_run_s {
	/* Index of bin this run is associated with. */
	szind_t		binind;

	/* Number of free regions in run. */
	unsigned	nfree;

	/* Per region allocated/deallocated bitmap. */
	bitmap_t	bitmap[BITMAP_GROUPS_MAX];
};

/* Each element of the chunk map corresponds to one page within the chunk. */
struct arena_chunk_map_bits_s {
	/*
	 * Run address (or size) and various flags are stored together.  The bit
	 * layout looks like (assuming 32-bit system):
	 *
	 *   ???????? ???????? ???nnnnn nnndumla
	 *
	 * ? : Unallocated: Run address for first/last pages, unset for internal
	 *                  pages.
	 *     Small: Run page offset.
	 *     Large: Run page count for first page, unset for trailing pages.
	 * n : binind for small size class, BININD_INVALID for large size class.
	 * d : dirty?
	 * u : unzeroed?
	 * m : decommitted?
	 * l : large?
	 * a : allocated?
	 *
	 * Following are example bit patterns for the three types of runs.
	 *
	 * p : run page offset
	 * s : run size
	 * n : binind for size class; large objects set these to BININD_INVALID
	 * x : don't care
	 * - : 0
	 * + : 1
	 * [DUMLA] : bit set
	 * [dumla] : bit unset
	 *
	 *   Unallocated (clean):
	 *     ssssssss ssssssss sss+++++ +++dum-a
	 *     xxxxxxxx xxxxxxxx xxxxxxxx xxx-Uxxx
	 *     ssssssss ssssssss sss+++++ +++dUm-a
	 *
	 *   Unallocated (dirty):
	 *     ssssssss ssssssss sss+++++ +++D-m-a
	 *     xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
	 *     ssssssss ssssssss sss+++++ +++D-m-a
	 *
	 *   Small:
	 *     pppppppp pppppppp pppnnnnn nnnd---A
	 *     pppppppp pppppppp pppnnnnn nnn----A
	 *     pppppppp pppppppp pppnnnnn nnnd---A
	 *
	 *   Large:
	 *     ssssssss ssssssss sss+++++ +++D--LA
	 *     xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
	 *     -------- -------- ---+++++ +++D--LA
	 *
	 *   Large (sampled, size <= LARGE_MINCLASS):
	 *     ssssssss ssssssss sssnnnnn nnnD--LA
	 *     xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
	 *     -------- -------- ---+++++ +++D--LA
	 *
	 *   Large (not sampled, size == LARGE_MINCLASS):
	 *     ssssssss ssssssss sss+++++ +++D--LA
	 *     xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
	 *     -------- -------- ---+++++ +++D--LA
	 */
	size_t				bits;
#define	CHUNK_MAP_ALLOCATED	((size_t)0x01U)
#define	CHUNK_MAP_LARGE		((size_t)0x02U)
#define	CHUNK_MAP_STATE_MASK	((size_t)0x3U)

#define	CHUNK_MAP_DECOMMITTED	((size_t)0x04U)
#define	CHUNK_MAP_UNZEROED	((size_t)0x08U)
#define	CHUNK_MAP_DIRTY		((size_t)0x10U)
#define	CHUNK_MAP_FLAGS_MASK	((size_t)0x1cU)

#define	CHUNK_MAP_BININD_SHIFT	5
#define	BININD_INVALID		((size_t)0xffU)
#define	CHUNK_MAP_BININD_MASK	(BININD_INVALID << CHUNK_MAP_BININD_SHIFT)
#define	CHUNK_MAP_BININD_INVALID CHUNK_MAP_BININD_MASK

#define	CHUNK_MAP_RUNIND_SHIFT	(CHUNK_MAP_BININD_SHIFT + 8)
#define	CHUNK_MAP_SIZE_SHIFT	(CHUNK_MAP_RUNIND_SHIFT - LG_PAGE)
#define	CHUNK_MAP_SIZE_MASK						\
    (~(CHUNK_MAP_BININD_MASK | CHUNK_MAP_FLAGS_MASK | CHUNK_MAP_STATE_MASK))
};

struct arena_runs_dirty_link_s {
	qr(arena_runs_dirty_link_t)	rd_link;
};

/*
 * Each arena_chunk_map_misc_t corresponds to one page within the chunk, just
 * like arena_chunk_map_bits_t.  Two separate arrays are stored within each
 * chunk header in order to improve cache locality.
 */
struct arena_chunk_map_misc_s {
	/*
	 * Linkage for run heaps.  There are two disjoint uses:
	 *
	 * 1) arena_t's runs_avail heaps.
	 * 2) arena_run_t conceptually uses this linkage for in-use non-full
	 *    runs, rather than directly embedding linkage.
	 */
	phn(arena_chunk_map_misc_t)		ph_link;

	union {
		/* Linkage for list of dirty runs. */
		arena_runs_dirty_link_t		rd;

		/* Profile counters, used for large object runs. */
		union {
			void			*prof_tctx_pun;
			prof_tctx_t		*prof_tctx;
		};

		/* Small region run metadata. */
		arena_run_t			run;
	};
};
typedef ph(arena_chunk_map_misc_t) arena_run_heap_t;
#endif /* JEMALLOC_ARENA_STRUCTS_A */

#ifdef JEMALLOC_ARENA_STRUCTS_B
/* Arena chunk header. */
struct arena_chunk_s {
	/*
	 * A pointer to the arena that owns the chunk is stored within the node.
	 * This field as a whole is used by chunks_rtree to support both
	 * ivsalloc() and core-based debugging.
	 */
	extent_node_t		node;

	/*
	 * Map of pages within chunk that keeps track of free/large/small.  The
	 * first map_bias entries are omitted, since the chunk header does not
	 * need to be tracked in the map.  This omission saves a header page
	 * for common chunk sizes (e.g. 4 MiB).
	 */
	arena_chunk_map_bits_t	map_bits[1]; /* Dynamically sized. */
};

/*
 * Read-only information associated with each element of arena_t's bins array
 * is stored separately, partly to reduce memory usage (only one copy, rather
 * than one per arena), but mainly to avoid false cacheline sharing.
 *
 * Each run has the following layout:
 *
 *               /--------------------\
 *               | pad?               |
 *               |--------------------|
 *               | redzone            |
 *   reg0_offset | region 0           |
 *               | redzone            |
 *               |--------------------| \
 *               | redzone            | |
 *               | region 1           |  > reg_interval
 *               | redzone            | /
 *               |--------------------|
 *               | ...                |
 *               | ...                |
 *               | ...                |
 *               |--------------------|
 *               | redzone            |
 *               | region nregs-1     |
 *               | redzone            |
 *               |--------------------|
 *               | alignment pad?     |
 *               \--------------------/
 *
 * reg_interval has at least the same minimum alignment as reg_size; this
 * preserves the alignment constraint that sa2u() depends on.  Alignment pad is
 * either 0 or redzone_size; it is present only if needed to align reg0_offset.
 */
struct arena_bin_info_s {
	/* Size of regions in a run for this bin's size class. */
	size_t			reg_size;

	/* Redzone size. */
	size_t			redzone_size;

	/* Interval between regions (reg_size + (redzone_size << 1)). */
	size_t			reg_interval;

	/* Total size of a run for this bin's size class. */
	size_t			run_size;

	/* Total number of regions in a run for this bin's size class. */
	uint32_t		nregs;

	/*
	 * Metadata used to manipulate bitmaps for runs associated with this
	 * bin.
	 */
	bitmap_info_t		bitmap_info;

	/* Offset of first region in a run for this bin's size class. */
	uint32_t		reg0_offset;
};

struct arena_bin_s {
	/*
	 * All operations on runcur, runs, and stats require that lock be
	 * locked.  Run allocation/deallocation are protected by the arena lock,
	 * which may be acquired while holding one or more bin locks, but not
	 * vise versa.
	 */
	malloc_mutex_t		lock;

	/*
	 * Current run being used to service allocations of this bin's size
	 * class.
	 */
	arena_run_t		*runcur;

	/*
	 * Heap of non-full runs.  This heap is used when looking for an
	 * existing run when runcur is no longer usable.  We choose the
	 * non-full run that is lowest in memory; this policy tends to keep
	 * objects packed well, and it can also help reduce the number of
	 * almost-empty chunks.
	 */
	arena_run_heap_t	runs;

	/* Bin statistics. */
	malloc_bin_stats_t	stats;
};

#define ACACHE_SIZE_RATIO_DEFAULT 4 /* default acache size = (4 * tcache) */
#define ACACHE_BYPASS_IND_DEFAULT 16

#define ACACHE_EPOCH_OFF 32
#define ACACHE_EPOCH_INC ((uint64_t)1 << ACACHE_EPOCH_OFF)
#define ACACHE_LOCKBIT   ((uint64_t)1 << 31)
#define ACACHE_NCACHED_BITS 15
#define ACACHE_NCACHED_MASK (((uint64_t)1 << ACACHE_NCACHED_BITS) - 1)
#define ACACHE_CAS_FAIL_RETRY_MAX 8

#define SPINWAIT_RETRY { CPU_SPINWAIT; goto label_retry; }

struct arena_cache_bin_info_s {
	unsigned	ncached_max;	/* Upper limit on ncached. */
	unsigned	flush_remain;
};

struct arena_cache_bin_s {
	/*
	 * +--------     Layout for data (64 bits)     --------+
	 * | 63...32 |  31  |   30   | 29.......15 | 14......0 |
	 * | [epoch] | lock | unused | [low_water] | [ncached] |
	 * +---------------------------------------------------+
	 * The single bit lock is used for cache_free.
	 */
	acache_state_t		data;
	void			**avail;	/* Stack of available objects. */

	/* Queued stats updated w/ atomics. */
	uint64_t		queued_nflushes;
	uint64_t		queued_nrequests;

	/*
	 * For small items in acache, we maintain the following info and only reuse
	 * these memory if it's freed by the same thread.
	 */
	tsdn_t			*last_thd; 				/* Most recent thread dalloced to acache */
	uint64_t		n_items_last_thd;	/* # of items freed by last_thd */
};

struct arena_cache_s {
	unsigned next_gc_bin;
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
	 * When perCPU arena is enabled, to amortize the cost of reading / updating
	 * the current CPU id, track the most recent thread accessing this arena, and
	 * only read CPU if there is a mismatch.
	 */
	tsdn_t		*last_thd;

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
	 * 3) Chunk- and run-related operations are protected by this mutex.
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


	/* Extant arena chunks. */
	ql_head(extent_node_t)	achunks;

	/*
	 * In order to avoid rapid chunk allocation/deallocation when an arena
	 * oscillates right on the cusp of needing a new chunk, cache the most
	 * recently freed chunk.  The spare is left in the arena's chunk trees
	 * until it is deleted.
	 *
	 * There is one spare chunk per arena, rather than one spare total, in
	 * order to avoid interactions between multiple threads that could make
	 * a single spare inadequate.
	 */
	arena_chunk_t		*spare;

	/* Minimum ratio (log base 2) of nactive:ndirty. */
	ssize_t			lg_dirty_mult;

	/* True if a thread is currently executing arena_purge_to_limit(). */
	bool			purging;

	/* Number of pages in active runs and huge regions. */
	size_t			nactive;

	/*
	 * Current count of pages within unused runs that are potentially
	 * dirty, and for which madvise(... MADV_DONTNEED) has not been called.
	 * By tracking this, we can institute a limit on how much dirty unused
	 * memory is mapped for each arena.
	 */
	size_t			ndirty;

	/*
	 * Unused dirty memory this arena manages.  Dirty memory is conceptually
	 * tracked as an arbitrarily interleaved LRU of dirty runs and cached
	 * chunks, but the list linkage is actually semi-duplicated in order to
	 * avoid extra arena_chunk_map_misc_t space overhead.
	 *
	 *   LRU-----------------------------------------------------------MRU
	 *
	 *        /-- arena ---\
	 *        |            |
	 *        |            |
	 *        |------------|                             /- chunk -\
	 *   ...->|chunks_cache|<--------------------------->|  /----\ |<--...
	 *        |------------|                             |  |node| |
	 *        |            |                             |  |    | |
	 *        |            |    /- run -\    /- run -\   |  |    | |
	 *        |            |    |       |    |       |   |  |    | |
	 *        |            |    |       |    |       |   |  |    | |
	 *        |------------|    |-------|    |-------|   |  |----| |
	 *   ...->|runs_dirty  |<-->|rd     |<-->|rd     |<---->|rd  |<----...
	 *        |------------|    |-------|    |-------|   |  |----| |
	 *        |            |    |       |    |       |   |  |    | |
	 *        |            |    |       |    |       |   |  \----/ |
	 *        |            |    \-------/    \-------/   |         |
	 *        |            |                             |         |
	 *        |            |                             |         |
	 *        \------------/                             \---------/
	 */
	arena_runs_dirty_link_t	runs_dirty;
	extent_node_t		chunks_cache;

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

	/* Extant huge allocations. */
	ql_head(extent_node_t)	huge;
	/* Synchronizes all huge allocation/update/deallocation. */
	malloc_mutex_t		huge_mtx;

	/*
	 * Trees of chunks that were previously allocated (trees differ only in
	 * node ordering).  These are used when allocating chunks, in an attempt
	 * to re-use address space.  Depending on function, different tree
	 * orderings are needed, which is why there are two trees with the same
	 * contents.
	 */
	extent_tree_t		chunks_szad_cached;
	extent_tree_t		chunks_ad_cached;
	extent_tree_t		chunks_szad_retained;
	extent_tree_t		chunks_ad_retained;

	malloc_mutex_t		chunks_mtx;
	/* Cache of nodes that were allocated via base_alloc(). */
	ql_head(extent_node_t)	node_cache;
	malloc_mutex_t		node_cache_mtx;

	/* User-configurable chunk hook functions. */
	chunk_hooks_t		chunk_hooks;

	/* bins is used to store trees of free regions. */
	arena_bin_t		bins[NBINS];

	pthread_t		purge_thread;
	/* Bins in acache is dynamically sized like tcache bins. When we get rid of
	 * runs_avail (last member currently), acache can be the last and doesn't have
	 * to be a pointer. */
	arena_cache_t *acache;

	/*
	 * Quantized address-ordered heaps of this arena's available runs.  The
	 * heaps are used for first-best-fit run allocation.
	 */
	arena_run_heap_t	runs_avail[1]; /* Dynamically sized. */
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

extern bool	opt_acache;
extern unsigned opt_acache_size_ratio;
extern unsigned opt_acache_bypass;

extern unsigned opt_percpu_arena;
extern bool	opt_arena_purging_thread;

extern purge_mode_t	opt_purge;
extern const char	*purge_mode_names[];
extern ssize_t		opt_lg_dirty_mult;
extern ssize_t		opt_decay_time;

extern arena_bin_info_t	arena_bin_info[NBINS];
extern arena_cache_bin_info_t	*arena_cache_bin_info;
extern unsigned	nhbins;

extern size_t		map_bias; /* Number of arena chunk header pages. */
extern size_t		map_misc_offset;
extern size_t		arena_maxrun; /* Max run size for arenas. */
extern size_t		large_maxclass; /* Max large size class. */
extern size_t		run_quantize_max; /* Max run_quantize_*() input. */
extern unsigned		nlclasses; /* Number of large size classes. */
extern unsigned		nhclasses; /* Number of huge size classes. */

#ifdef JEMALLOC_JET
typedef size_t (run_quantize_t)(size_t);
extern run_quantize_t *run_quantize_floor;
extern run_quantize_t *run_quantize_ceil;
#endif
void	arena_chunk_cache_maybe_insert(arena_t *arena, extent_node_t *node,
    bool cache);
void	arena_chunk_cache_maybe_remove(arena_t *arena, extent_node_t *node,
    bool cache);
extent_node_t	*arena_node_alloc(tsdn_t *tsdn, arena_t *arena);
void	arena_node_dalloc(tsdn_t *tsdn, arena_t *arena, extent_node_t *node);
void	*arena_chunk_alloc_huge(tsdn_t *tsdn, arena_t *arena, size_t usize,
    size_t alignment, bool *zero);
void	arena_chunk_dalloc_huge(tsdn_t *tsdn, arena_t *arena, void *chunk,
    size_t usize);
void	arena_chunk_ralloc_huge_similar(tsdn_t *tsdn, arena_t *arena,
    void *chunk, size_t oldsize, size_t usize);
void	arena_chunk_ralloc_huge_shrink(tsdn_t *tsdn, arena_t *arena,
    void *chunk, size_t oldsize, size_t usize);
bool	arena_chunk_ralloc_huge_expand(tsdn_t *tsdn, arena_t *arena,
    void *chunk, size_t oldsize, size_t usize, bool *zero);
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
void	arena_alloc_junk_small(void *ptr, arena_bin_info_t *bin_info,
    bool zero);
#ifdef JEMALLOC_JET
typedef void (arena_redzone_corruption_t)(void *, size_t, bool, size_t,
    uint8_t);
extern arena_redzone_corruption_t *arena_redzone_corruption;
typedef void (arena_dalloc_junk_small_t)(void *, arena_bin_info_t *);
extern arena_dalloc_junk_small_t *arena_dalloc_junk_small;
#else
void	arena_dalloc_junk_small(void *ptr, arena_bin_info_t *bin_info);
#endif
void	arena_quarantine_junk_small(void *ptr, size_t usize);
void	*arena_malloc_large(tsdn_t *tsdn, arena_t *arena, szind_t ind,
    bool zero);
void	*arena_malloc_hard(tsdn_t *tsdn, arena_t *arena, size_t size,
    szind_t ind, bool zero);
void	*arena_palloc(tsdn_t *tsdn, arena_t *arena, size_t usize,
    size_t alignment, bool zero, tcache_t *tcache);
void	arena_prof_promoted(tsdn_t *tsdn, const void *ptr, size_t size);
void	arena_dalloc_bin_junked_locked(tsdn_t *tsdn, arena_t *arena,
    arena_chunk_t *chunk, void *ptr, arena_chunk_map_bits_t *bitselm);
void	arena_dalloc_bin(tsdn_t *tsdn, arena_t *arena, arena_chunk_t *chunk,
    void *ptr, size_t pageind, arena_chunk_map_bits_t *bitselm);
void	arena_dalloc_small(tsdn_t *tsdn, arena_t *arena, arena_chunk_t *chunk,
    void *ptr, size_t pageind);
#ifdef JEMALLOC_JET
typedef void (arena_dalloc_junk_large_t)(void *, size_t);
extern arena_dalloc_junk_large_t *arena_dalloc_junk_large;
#else
void	arena_dalloc_junk_large(void *ptr, size_t usize);
#endif
void	arena_dalloc_large_junked_locked(tsdn_t *tsdn, arena_t *arena,
    arena_chunk_t *chunk, void *ptr);
void	arena_dalloc_large(tsdn_t *tsdn, arena_t *arena, arena_chunk_t *chunk,
    void *ptr);
#ifdef JEMALLOC_JET
typedef void (arena_ralloc_junk_large_t)(void *, size_t, size_t);
extern arena_ralloc_junk_large_t *arena_ralloc_junk_large;
#endif
bool	arena_ralloc_no_move(tsdn_t *tsdn, void *ptr, size_t oldsize,
    size_t size, size_t extra, bool zero);
void	*arena_ralloc(tsd_t *tsd, arena_t *arena, void *ptr, size_t oldsize,
    size_t size, size_t alignment, bool zero, tcache_t *tcache);
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
    malloc_bin_stats_t *bstats, malloc_large_stats_t *lstats,
    malloc_huge_stats_t *hstats);
unsigned	arena_nthreads_get(arena_t *arena, bool internal);
void	arena_nthreads_inc(arena_t *arena, bool internal);
void	arena_nthreads_dec(arena_t *arena, bool internal);
arena_t	*arena_new(tsdn_t *tsdn, unsigned ind);
bool	arena_boot(void);
void	arena_prefork0(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork1(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork2(tsdn_t *tsdn, arena_t *arena);
void	arena_prefork3(tsdn_t *tsdn, arena_t *arena);
void	arena_postfork_parent(tsdn_t *tsdn, arena_t *arena);
void	arena_postfork_child(tsdn_t *tsdn, arena_t *arena);
bool	a0_purge_thread_init(void);


#endif /* JEMALLOC_H_EXTERNS */
/******************************************************************************/
#ifdef JEMALLOC_H_INLINES

#ifndef JEMALLOC_ENABLE_INLINE
arena_chunk_map_bits_t	*arena_bitselm_get_mutable(arena_chunk_t *chunk,
    size_t pageind);
const arena_chunk_map_bits_t	*arena_bitselm_get_const(
    const arena_chunk_t *chunk, size_t pageind);
arena_chunk_map_misc_t	*arena_miscelm_get_mutable(arena_chunk_t *chunk,
    size_t pageind);
const arena_chunk_map_misc_t	*arena_miscelm_get_const(
    const arena_chunk_t *chunk, size_t pageind);
size_t	arena_miscelm_to_pageind(const arena_chunk_map_misc_t *miscelm);
void	*arena_miscelm_to_rpages(const arena_chunk_map_misc_t *miscelm);
arena_chunk_map_misc_t	*arena_rd_to_miscelm(arena_runs_dirty_link_t *rd);
arena_chunk_map_misc_t	*arena_run_to_miscelm(arena_run_t *run);
size_t	*arena_mapbitsp_get_mutable(arena_chunk_t *chunk, size_t pageind);
const size_t	*arena_mapbitsp_get_const(const arena_chunk_t *chunk,
    size_t pageind);
size_t	arena_mapbitsp_read(const size_t *mapbitsp);
size_t	arena_mapbits_get(const arena_chunk_t *chunk, size_t pageind);
size_t	arena_mapbits_size_decode(size_t mapbits);
size_t	arena_mapbits_unallocated_size_get(const arena_chunk_t *chunk,
    size_t pageind);
size_t	arena_mapbits_large_size_get(const arena_chunk_t *chunk,
    size_t pageind);
size_t	arena_mapbits_small_runind_get(const arena_chunk_t *chunk,
    size_t pageind);
szind_t	arena_mapbits_binind_get(const arena_chunk_t *chunk, size_t pageind);
size_t	arena_mapbits_dirty_get(const arena_chunk_t *chunk, size_t pageind);
size_t	arena_mapbits_unzeroed_get(const arena_chunk_t *chunk, size_t pageind);
size_t	arena_mapbits_decommitted_get(const arena_chunk_t *chunk,
    size_t pageind);
size_t	arena_mapbits_large_get(const arena_chunk_t *chunk, size_t pageind);
size_t	arena_mapbits_allocated_get(const arena_chunk_t *chunk, size_t pageind);
void	arena_mapbitsp_write(size_t *mapbitsp, size_t mapbits);
size_t	arena_mapbits_size_encode(size_t size);
void	arena_mapbits_unallocated_set(arena_chunk_t *chunk, size_t pageind,
    size_t size, size_t flags);
void	arena_mapbits_unallocated_size_set(arena_chunk_t *chunk, size_t pageind,
    size_t size);
void	arena_mapbits_internal_set(arena_chunk_t *chunk, size_t pageind,
    size_t flags);
void	arena_mapbits_large_set(arena_chunk_t *chunk, size_t pageind,
    size_t size, size_t flags);
void	arena_mapbits_large_binind_set(arena_chunk_t *chunk, size_t pageind,
    szind_t binind);
void	arena_mapbits_small_set(arena_chunk_t *chunk, size_t pageind,
    size_t runind, szind_t binind, size_t flags);
void	arena_metadata_allocated_add(arena_t *arena, size_t size);
void	arena_metadata_allocated_sub(arena_t *arena, size_t size);
size_t	arena_metadata_allocated_get(arena_t *arena);
bool	arena_prof_accum_impl(arena_t *arena, uint64_t accumbytes);
bool	arena_prof_accum_locked(arena_t *arena, uint64_t accumbytes);
bool	arena_prof_accum(tsdn_t *tsdn, arena_t *arena, uint64_t accumbytes);
szind_t	arena_ptr_small_binind_get(const void *ptr, size_t mapbits);
szind_t	arena_bin_index(arena_t *arena, arena_bin_t *bin);
size_t	arena_run_regind(arena_run_t *run, arena_bin_info_t *bin_info,
    const void *ptr);
prof_tctx_t	*arena_prof_tctx_get(tsdn_t *tsdn, const void *ptr);
void	arena_prof_tctx_set(tsdn_t *tsdn, const void *ptr, size_t usize,
    prof_tctx_t *tctx);
void	arena_prof_tctx_reset(tsdn_t *tsdn, const void *ptr, size_t usize,
    const void *old_ptr, prof_tctx_t *old_tctx);
void	arena_decay_ticks(tsdn_t *tsdn, arena_t *arena, unsigned nticks);
void	arena_decay_tick(tsdn_t *tsdn, arena_t *arena);
void	*arena_malloc(tsdn_t *tsdn, arena_t *arena, size_t size, szind_t ind,
    bool zero, tcache_t *tcache, bool slow_path);
arena_t	*arena_aalloc(const void *ptr);
size_t	arena_salloc(tsdn_t *tsdn, const void *ptr, bool demote);
void	arena_dalloc(tsdn_t *tsdn, void *ptr, tcache_t *tcache, bool slow_path);
void	arena_sdalloc(tsdn_t *tsdn, void *ptr, size_t size, tcache_t *tcache,
    bool slow_path);
size_t arena_cache_alloc_small(tsdn_t *tsdn, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, szind_t binind);
void *arena_cache_alloc_large(tsdn_t *tsdn, arena_t *arena, szind_t binind,
    size_t usize, bool zero);
void arena_cache_dalloc(tsdn_t *tsdn, arena_t *arena, void **items,
    size_t n_items, szind_t binind, uint64_t nrequests, const bool is_large);
void arena_cache_flush(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
    arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
    uint64_t nrequests, const bool is_large);
void arena_cache_gc(tsdn_t *tsdn, arena_t *arena);
void arena_cache_merge_stats(arena_t *arena, szind_t binind,
    uint64_t nrequests);
#endif

#if (defined(JEMALLOC_ENABLE_INLINE) || defined(JEMALLOC_ARENA_C_))
#  ifdef JEMALLOC_ARENA_INLINE_A
JEMALLOC_ALWAYS_INLINE arena_chunk_map_bits_t *
arena_bitselm_get_mutable(arena_chunk_t *chunk, size_t pageind)
{

	assert(pageind >= map_bias);
	assert(pageind < chunk_npages);

	return (&chunk->map_bits[pageind-map_bias]);
}

JEMALLOC_ALWAYS_INLINE const arena_chunk_map_bits_t *
arena_bitselm_get_const(const arena_chunk_t *chunk, size_t pageind)
{

	return (arena_bitselm_get_mutable((arena_chunk_t *)chunk, pageind));
}

JEMALLOC_ALWAYS_INLINE arena_chunk_map_misc_t *
arena_miscelm_get_mutable(arena_chunk_t *chunk, size_t pageind)
{

	assert(pageind >= map_bias);
	assert(pageind < chunk_npages);

	return ((arena_chunk_map_misc_t *)((uintptr_t)chunk +
	    (uintptr_t)map_misc_offset) + pageind-map_bias);
}

JEMALLOC_ALWAYS_INLINE const arena_chunk_map_misc_t *
arena_miscelm_get_const(const arena_chunk_t *chunk, size_t pageind)
{

	return (arena_miscelm_get_mutable((arena_chunk_t *)chunk, pageind));
}

JEMALLOC_ALWAYS_INLINE size_t
arena_miscelm_to_pageind(const arena_chunk_map_misc_t *miscelm)
{
	arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(miscelm);
	size_t pageind = ((uintptr_t)miscelm - ((uintptr_t)chunk +
	    map_misc_offset)) / sizeof(arena_chunk_map_misc_t) + map_bias;

	assert(pageind >= map_bias);
	assert(pageind < chunk_npages);

	return (pageind);
}

JEMALLOC_ALWAYS_INLINE void *
arena_miscelm_to_rpages(const arena_chunk_map_misc_t *miscelm)
{
	arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(miscelm);
	size_t pageind = arena_miscelm_to_pageind(miscelm);

	return ((void *)((uintptr_t)chunk + (pageind << LG_PAGE)));
}

JEMALLOC_ALWAYS_INLINE arena_chunk_map_misc_t *
arena_rd_to_miscelm(arena_runs_dirty_link_t *rd)
{
	arena_chunk_map_misc_t *miscelm = (arena_chunk_map_misc_t
	    *)((uintptr_t)rd - offsetof(arena_chunk_map_misc_t, rd));

	assert(arena_miscelm_to_pageind(miscelm) >= map_bias);
	assert(arena_miscelm_to_pageind(miscelm) < chunk_npages);

	return (miscelm);
}

JEMALLOC_ALWAYS_INLINE arena_chunk_map_misc_t *
arena_run_to_miscelm(arena_run_t *run)
{
	arena_chunk_map_misc_t *miscelm = (arena_chunk_map_misc_t
	    *)((uintptr_t)run - offsetof(arena_chunk_map_misc_t, run));

	assert(arena_miscelm_to_pageind(miscelm) >= map_bias);
	assert(arena_miscelm_to_pageind(miscelm) < chunk_npages);

	return (miscelm);
}

JEMALLOC_ALWAYS_INLINE size_t *
arena_mapbitsp_get_mutable(arena_chunk_t *chunk, size_t pageind)
{

	return (&arena_bitselm_get_mutable(chunk, pageind)->bits);
}

JEMALLOC_ALWAYS_INLINE const size_t *
arena_mapbitsp_get_const(const arena_chunk_t *chunk, size_t pageind)
{

	return (arena_mapbitsp_get_mutable((arena_chunk_t *)chunk, pageind));
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbitsp_read(const size_t *mapbitsp)
{

	return (*mapbitsp);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_get(const arena_chunk_t *chunk, size_t pageind)
{

	return (arena_mapbitsp_read(arena_mapbitsp_get_const(chunk, pageind)));
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_size_decode(size_t mapbits)
{
	size_t size;

#if CHUNK_MAP_SIZE_SHIFT > 0
	size = (mapbits & CHUNK_MAP_SIZE_MASK) >> CHUNK_MAP_SIZE_SHIFT;
#elif CHUNK_MAP_SIZE_SHIFT == 0
	size = mapbits & CHUNK_MAP_SIZE_MASK;
#else
	size = (mapbits & CHUNK_MAP_SIZE_MASK) << -CHUNK_MAP_SIZE_SHIFT;
#endif

	return (size);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_unallocated_size_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & (CHUNK_MAP_LARGE|CHUNK_MAP_ALLOCATED)) == 0);
	return (arena_mapbits_size_decode(mapbits));
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_large_size_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & (CHUNK_MAP_LARGE|CHUNK_MAP_ALLOCATED)) ==
	    (CHUNK_MAP_LARGE|CHUNK_MAP_ALLOCATED));
	return (arena_mapbits_size_decode(mapbits));
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_small_runind_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & (CHUNK_MAP_LARGE|CHUNK_MAP_ALLOCATED)) ==
	    CHUNK_MAP_ALLOCATED);
	return (mapbits >> CHUNK_MAP_RUNIND_SHIFT);
}

JEMALLOC_ALWAYS_INLINE szind_t
arena_mapbits_binind_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;
	szind_t binind;

	mapbits = arena_mapbits_get(chunk, pageind);
	binind = (mapbits & CHUNK_MAP_BININD_MASK) >> CHUNK_MAP_BININD_SHIFT;
	assert(binind < NBINS || binind == BININD_INVALID);
	return (binind);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_dirty_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & CHUNK_MAP_DECOMMITTED) == 0 || (mapbits &
	    (CHUNK_MAP_DIRTY|CHUNK_MAP_UNZEROED)) == 0);
	return (mapbits & CHUNK_MAP_DIRTY);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_unzeroed_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & CHUNK_MAP_DECOMMITTED) == 0 || (mapbits &
	    (CHUNK_MAP_DIRTY|CHUNK_MAP_UNZEROED)) == 0);
	return (mapbits & CHUNK_MAP_UNZEROED);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_decommitted_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	assert((mapbits & CHUNK_MAP_DECOMMITTED) == 0 || (mapbits &
	    (CHUNK_MAP_DIRTY|CHUNK_MAP_UNZEROED)) == 0);
	return (mapbits & CHUNK_MAP_DECOMMITTED);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_large_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	return (mapbits & CHUNK_MAP_LARGE);
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_allocated_get(const arena_chunk_t *chunk, size_t pageind)
{
	size_t mapbits;

	mapbits = arena_mapbits_get(chunk, pageind);
	return (mapbits & CHUNK_MAP_ALLOCATED);
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbitsp_write(size_t *mapbitsp, size_t mapbits)
{

	*mapbitsp = mapbits;
}

JEMALLOC_ALWAYS_INLINE size_t
arena_mapbits_size_encode(size_t size)
{
	size_t mapbits;

#if CHUNK_MAP_SIZE_SHIFT > 0
	mapbits = size << CHUNK_MAP_SIZE_SHIFT;
#elif CHUNK_MAP_SIZE_SHIFT == 0
	mapbits = size;
#else
	mapbits = size >> -CHUNK_MAP_SIZE_SHIFT;
#endif

	assert((mapbits & ~CHUNK_MAP_SIZE_MASK) == 0);
	return (mapbits);
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_unallocated_set(arena_chunk_t *chunk, size_t pageind, size_t size,
    size_t flags)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);

	assert((size & PAGE_MASK) == 0);
	assert((flags & CHUNK_MAP_FLAGS_MASK) == flags);
	assert((flags & CHUNK_MAP_DECOMMITTED) == 0 || (flags &
	    (CHUNK_MAP_DIRTY|CHUNK_MAP_UNZEROED)) == 0);
	arena_mapbitsp_write(mapbitsp, arena_mapbits_size_encode(size) |
	    CHUNK_MAP_BININD_INVALID | flags);
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_unallocated_size_set(arena_chunk_t *chunk, size_t pageind,
    size_t size)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);
	size_t mapbits = arena_mapbitsp_read(mapbitsp);

	assert((size & PAGE_MASK) == 0);
	assert((mapbits & (CHUNK_MAP_LARGE|CHUNK_MAP_ALLOCATED)) == 0);
	arena_mapbitsp_write(mapbitsp, arena_mapbits_size_encode(size) |
	    (mapbits & ~CHUNK_MAP_SIZE_MASK));
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_internal_set(arena_chunk_t *chunk, size_t pageind, size_t flags)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);

	assert((flags & CHUNK_MAP_UNZEROED) == flags);
	arena_mapbitsp_write(mapbitsp, flags);
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_large_set(arena_chunk_t *chunk, size_t pageind, size_t size,
    size_t flags)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);

	assert((size & PAGE_MASK) == 0);
	assert((flags & CHUNK_MAP_FLAGS_MASK) == flags);
	assert((flags & CHUNK_MAP_DECOMMITTED) == 0 || (flags &
	    (CHUNK_MAP_DIRTY|CHUNK_MAP_UNZEROED)) == 0);
	arena_mapbitsp_write(mapbitsp, arena_mapbits_size_encode(size) |
	    CHUNK_MAP_BININD_INVALID | flags | CHUNK_MAP_LARGE |
	    CHUNK_MAP_ALLOCATED);
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_large_binind_set(arena_chunk_t *chunk, size_t pageind,
    szind_t binind)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);
	size_t mapbits = arena_mapbitsp_read(mapbitsp);

	assert(binind <= BININD_INVALID);
	assert(arena_mapbits_large_size_get(chunk, pageind) == LARGE_MINCLASS +
	    large_pad);
	arena_mapbitsp_write(mapbitsp, (mapbits & ~CHUNK_MAP_BININD_MASK) |
	    (binind << CHUNK_MAP_BININD_SHIFT));
}

JEMALLOC_ALWAYS_INLINE void
arena_mapbits_small_set(arena_chunk_t *chunk, size_t pageind, size_t runind,
    szind_t binind, size_t flags)
{
	size_t *mapbitsp = arena_mapbitsp_get_mutable(chunk, pageind);

	assert(binind < BININD_INVALID);
	assert(pageind - runind >= map_bias);
	assert((flags & CHUNK_MAP_UNZEROED) == flags);
	arena_mapbitsp_write(mapbitsp, (runind << CHUNK_MAP_RUNIND_SHIFT) |
	    (binind << CHUNK_MAP_BININD_SHIFT) | flags | CHUNK_MAP_ALLOCATED);
}

JEMALLOC_INLINE void
arena_metadata_allocated_add(arena_t *arena, size_t size)
{

	atomic_add_z(&arena->stats.metadata_allocated, size);
}

JEMALLOC_INLINE void
arena_metadata_allocated_sub(arena_t *arena, size_t size)
{

	atomic_sub_z(&arena->stats.metadata_allocated, size);
}

JEMALLOC_INLINE size_t
arena_metadata_allocated_get(arena_t *arena)
{

	return (atomic_read_z(&arena->stats.metadata_allocated));
}

JEMALLOC_INLINE bool
arena_prof_accum_impl(arena_t *arena, uint64_t accumbytes)
{

	cassert(config_prof);
	assert(prof_interval != 0);

	arena->prof_accumbytes += accumbytes;
	if (arena->prof_accumbytes >= prof_interval) {
		arena->prof_accumbytes -= prof_interval;
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
		uint64_t bytes_updated;

		bytes_updated = atomic_add_uint64(&arena->prof_accumbytes, accumbytes);

		if (bytes_updated >= prof_interval
		    && (bytes_updated - accumbytes < prof_interval)) {
			atomic_sub_uint64(&arena->prof_accumbytes, prof_interval);
			return (true);
		}
		return (false);
	}
}

JEMALLOC_ALWAYS_INLINE szind_t
arena_ptr_small_binind_get(const void *ptr, size_t mapbits)
{
	szind_t binind;

	binind = (mapbits & CHUNK_MAP_BININD_MASK) >> CHUNK_MAP_BININD_SHIFT;

	if (config_debug) {
		arena_chunk_t *chunk;
		arena_t *arena;
		size_t pageind;
		size_t actual_mapbits;
		size_t rpages_ind;
		const arena_run_t *run;
		arena_bin_t *bin;
		szind_t run_binind, actual_binind;
		arena_bin_info_t *bin_info;
		const arena_chunk_map_misc_t *miscelm;
		const void *rpages;

		assert(binind != BININD_INVALID);
		assert(binind < NBINS);
		chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
		arena = extent_node_arena_get(&chunk->node);
		pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >> LG_PAGE;
		actual_mapbits = arena_mapbits_get(chunk, pageind);
		assert(mapbits == actual_mapbits);
		assert(arena_mapbits_large_get(chunk, pageind) == 0);
		assert(arena_mapbits_allocated_get(chunk, pageind) != 0);
		rpages_ind = pageind - arena_mapbits_small_runind_get(chunk,
		    pageind);
		miscelm = arena_miscelm_get_const(chunk, rpages_ind);
		run = &miscelm->run;
		run_binind = run->binind;
		bin = &arena->bins[run_binind];
		actual_binind = (szind_t)(bin - arena->bins);
		assert(run_binind == actual_binind);
		bin_info = &arena_bin_info[actual_binind];
		rpages = arena_miscelm_to_rpages(miscelm);
		assert(((uintptr_t)ptr - ((uintptr_t)rpages +
		    (uintptr_t)bin_info->reg0_offset)) % bin_info->reg_interval
		    == 0);
	}

	return (binind);
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

JEMALLOC_INLINE size_t
arena_run_regind(arena_run_t *run, arena_bin_info_t *bin_info, const void *ptr)
{
	size_t diff, interval, shift, regind;
	arena_chunk_map_misc_t *miscelm = arena_run_to_miscelm(run);
	void *rpages = arena_miscelm_to_rpages(miscelm);

	/*
	 * Freeing a pointer lower than region zero can cause assertion
	 * failure.
	 */
	assert((uintptr_t)ptr >= (uintptr_t)rpages +
	    (uintptr_t)bin_info->reg0_offset);

	/*
	 * Avoid doing division with a variable divisor if possible.  Using
	 * actual division here can reduce allocator throughput by over 20%!
	 */
	diff = (size_t)((uintptr_t)ptr - (uintptr_t)rpages -
	    bin_info->reg0_offset);

	/* Rescale (factor powers of 2 out of the numerator and denominator). */
	interval = bin_info->reg_interval;
	shift = ffs_zu(interval) - 1;
	diff >>= shift;
	interval >>= shift;

	if (interval == 1) {
		/* The divisor was a power of 2. */
		regind = diff;
	} else {
		/*
		 * To divide by a number D that is not a power of two we
		 * multiply by (2^21 / D) and then right shift by 21 positions.
		 *
		 *   X / D
		 *
		 * becomes
		 *
		 *   (X * interval_invs[D - 3]) >> SIZE_INV_SHIFT
		 *
		 * We can omit the first three elements, because we never
		 * divide by 0, and 1 and 2 are both powers of two, which are
		 * handled above.
		 */
#define	SIZE_INV_SHIFT	((sizeof(size_t) << 3) - LG_RUN_MAXREGS)
#define	SIZE_INV(s)	(((ZU(1) << SIZE_INV_SHIFT) / (s)) + 1)
		static const size_t interval_invs[] = {
		    SIZE_INV(3),
		    SIZE_INV(4), SIZE_INV(5), SIZE_INV(6), SIZE_INV(7),
		    SIZE_INV(8), SIZE_INV(9), SIZE_INV(10), SIZE_INV(11),
		    SIZE_INV(12), SIZE_INV(13), SIZE_INV(14), SIZE_INV(15),
		    SIZE_INV(16), SIZE_INV(17), SIZE_INV(18), SIZE_INV(19),
		    SIZE_INV(20), SIZE_INV(21), SIZE_INV(22), SIZE_INV(23),
		    SIZE_INV(24), SIZE_INV(25), SIZE_INV(26), SIZE_INV(27),
		    SIZE_INV(28), SIZE_INV(29), SIZE_INV(30), SIZE_INV(31)
		};

		if (likely(interval <= ((sizeof(interval_invs) / sizeof(size_t))
		    + 2))) {
			regind = (diff * interval_invs[interval - 3]) >>
			    SIZE_INV_SHIFT;
		} else
			regind = diff / interval;
#undef SIZE_INV
#undef SIZE_INV_SHIFT
	}
	assert(diff == regind * interval);
	assert(regind < bin_info->nregs);

	return (regind);
}

JEMALLOC_INLINE prof_tctx_t *
arena_prof_tctx_get(tsdn_t *tsdn, const void *ptr)
{
	prof_tctx_t *ret;
	arena_chunk_t *chunk;

	cassert(config_prof);
	assert(ptr != NULL);

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr)) {
		size_t pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >> LG_PAGE;
		size_t mapbits = arena_mapbits_get(chunk, pageind);
		assert((mapbits & CHUNK_MAP_ALLOCATED) != 0);
		if (likely((mapbits & CHUNK_MAP_LARGE) == 0))
			ret = (prof_tctx_t *)(uintptr_t)1U;
		else {
			arena_chunk_map_misc_t *elm =
			    arena_miscelm_get_mutable(chunk, pageind);
			ret = atomic_read_p(&elm->prof_tctx_pun);
		}
	} else
		ret = huge_prof_tctx_get(tsdn, ptr);

	return (ret);
}

JEMALLOC_INLINE void
arena_prof_tctx_set(tsdn_t *tsdn, const void *ptr, size_t usize,
    prof_tctx_t *tctx)
{
	arena_chunk_t *chunk;

	cassert(config_prof);
	assert(ptr != NULL);

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr)) {
		size_t pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >> LG_PAGE;

		assert(arena_mapbits_allocated_get(chunk, pageind) != 0);

		if (unlikely(usize > SMALL_MAXCLASS || (uintptr_t)tctx >
		    (uintptr_t)1U)) {
			arena_chunk_map_misc_t *elm;

			assert(arena_mapbits_large_get(chunk, pageind) != 0);

			elm = arena_miscelm_get_mutable(chunk, pageind);
			atomic_write_p(&elm->prof_tctx_pun, tctx);
		} else {
			/*
			 * tctx must always be initialized for large runs.
			 * Assert that the surrounding conditional logic is
			 * equivalent to checking whether ptr refers to a large
			 * run.
			 */
			assert(arena_mapbits_large_get(chunk, pageind) == 0);
		}
	} else
		huge_prof_tctx_set(tsdn, ptr, tctx);
}

JEMALLOC_INLINE void
arena_prof_tctx_reset(tsdn_t *tsdn, const void *ptr, size_t usize,
    const void *old_ptr, prof_tctx_t *old_tctx)
{

	cassert(config_prof);
	assert(ptr != NULL);

	if (unlikely(usize > SMALL_MAXCLASS || (ptr == old_ptr &&
	    (uintptr_t)old_tctx > (uintptr_t)1U))) {
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
		if (likely(chunk != ptr)) {
			size_t pageind;
			arena_chunk_map_misc_t *elm;

			pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >>
			    LG_PAGE;
			assert(arena_mapbits_allocated_get(chunk, pageind) !=
			    0);
			assert(arena_mapbits_large_get(chunk, pageind) != 0);

			elm = arena_miscelm_get_mutable(chunk, pageind);
			atomic_write_p(&elm->prof_tctx_pun,
			    (prof_tctx_t *)(uintptr_t)1U);
		} else
			huge_prof_tctx_reset(tsdn, ptr);
	}
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

#ifdef JEMALLOC_HAVE_PTHREAD
	if (opt_arena_purging_thread)
		return;
#endif
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
arena_aalloc(const void *ptr)
{
	arena_chunk_t *chunk;

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr))
		return (extent_node_arena_get(&chunk->node));
	else
		return (huge_aalloc(ptr));
}

/* Return the size of the allocation pointed to by ptr. */
JEMALLOC_ALWAYS_INLINE size_t
arena_salloc(tsdn_t *tsdn, const void *ptr, bool demote)
{
	size_t ret;
	arena_chunk_t *chunk;
	size_t pageind;
	szind_t binind;

	assert(ptr != NULL);

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr)) {
		pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >> LG_PAGE;
		assert(arena_mapbits_allocated_get(chunk, pageind) != 0);
		binind = arena_mapbits_binind_get(chunk, pageind);
		if (unlikely(binind == BININD_INVALID || (config_prof && !demote
		    && arena_mapbits_large_get(chunk, pageind) != 0))) {
			/*
			 * Large allocation.  In the common case (demote), and
			 * as this is an inline function, most callers will only
			 * end up looking at binind to determine that ptr is a
			 * small allocation.
			 */
			assert(config_cache_oblivious || ((uintptr_t)ptr &
			    PAGE_MASK) == 0);
			ret = arena_mapbits_large_size_get(chunk, pageind) -
			    large_pad;
			assert(ret != 0);
			assert(pageind + ((ret+large_pad)>>LG_PAGE) <=
			    chunk_npages);
			assert(arena_mapbits_dirty_get(chunk, pageind) ==
			    arena_mapbits_dirty_get(chunk,
			    pageind+((ret+large_pad)>>LG_PAGE)-1));
		} else {
			/*
			 * Small allocation (possibly promoted to a large
			 * object).
			 */
			assert(arena_mapbits_large_get(chunk, pageind) != 0 ||
			    arena_ptr_small_binind_get(ptr,
			    arena_mapbits_get(chunk, pageind)) == binind);
			ret = index2size(binind);
		}
	} else
		ret = huge_salloc(tsdn, ptr);

	return (ret);
}

JEMALLOC_ALWAYS_INLINE void
arena_dalloc(tsdn_t *tsdn, void *ptr, tcache_t *tcache, bool slow_path)
{
	arena_chunk_t *chunk;
	size_t pageind, mapbits;

	assert(!tsdn_null(tsdn) || tcache == NULL);
	assert(ptr != NULL);

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr)) {
		pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >> LG_PAGE;
		mapbits = arena_mapbits_get(chunk, pageind);
		assert(arena_mapbits_allocated_get(chunk, pageind) != 0);
		if (likely((mapbits & CHUNK_MAP_LARGE) == 0)) {
			/* Small allocation. */
			if (likely(tcache != NULL)) {
				szind_t binind = arena_ptr_small_binind_get(ptr,
				    mapbits);
				tcache_dalloc_small(tsdn_tsd(tsdn), tcache, ptr,
				    binind, slow_path);
			} else {
				arena_dalloc_small(tsdn,
				    extent_node_arena_get(&chunk->node), chunk,
				    ptr, pageind);
			}
		} else {
			size_t size = arena_mapbits_large_size_get(chunk,
			    pageind);

			assert(config_cache_oblivious || ((uintptr_t)ptr &
			    PAGE_MASK) == 0);

			if (likely(tcache != NULL) && size - large_pad <=
			    tcache_maxclass) {
				tcache_dalloc_large(tsdn_tsd(tsdn), tcache, ptr,
				    size - large_pad, slow_path);
			} else {
				arena_dalloc_large(tsdn,
				    extent_node_arena_get(&chunk->node), chunk,
				    ptr);
			}
		}
	} else
		huge_dalloc(tsdn, ptr);
}

JEMALLOC_ALWAYS_INLINE void
arena_sdalloc(tsdn_t *tsdn, void *ptr, size_t size, tcache_t *tcache,
    bool slow_path)
{
	arena_chunk_t *chunk;

	assert(!tsdn_null(tsdn) || tcache == NULL);

	chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
	if (likely(chunk != ptr)) {
		if (config_prof && opt_prof) {
			size_t pageind = ((uintptr_t)ptr - (uintptr_t)chunk) >>
			    LG_PAGE;
			assert(arena_mapbits_allocated_get(chunk, pageind) !=
			    0);
			if (arena_mapbits_large_get(chunk, pageind) != 0) {
				/*
				 * Make sure to use promoted size, not request
				 * size.
				 */
				size = arena_mapbits_large_size_get(chunk,
				    pageind) - large_pad;
			}
		}
		assert(s2u(size) == s2u(arena_salloc(tsdn, ptr, false)));

		if (likely(size <= SMALL_MAXCLASS)) {
			/* Small allocation. */
			if (likely(tcache != NULL)) {
				szind_t binind = size2index(size);
				tcache_dalloc_small(tsdn_tsd(tsdn), tcache, ptr,
				    binind, slow_path);
			} else {
				size_t pageind = ((uintptr_t)ptr -
				    (uintptr_t)chunk) >> LG_PAGE;
				arena_dalloc_small(tsdn,
				    extent_node_arena_get(&chunk->node), chunk,
				    ptr, pageind);
			}
		} else {
			assert(config_cache_oblivious || ((uintptr_t)ptr &
			    PAGE_MASK) == 0);

			if (likely(tcache != NULL) && size <= tcache_maxclass) {
				tcache_dalloc_large(tsdn_tsd(tsdn), tcache, ptr,
				    size, slow_path);
			} else {
				arena_dalloc_large(tsdn,
				    extent_node_arena_get(&chunk->node), chunk,
				    ptr);
			}
		}
	} else
		huge_dalloc(tsdn, ptr);
}

JEMALLOC_ALWAYS_INLINE acache_state_t
cbin_state_get(arena_cache_bin_t *cbin, size_t *ncached, size_t *low_water,
    bool *locked)
{
	acache_state_t state = ACCESS_ONCE(cbin->data);

	if (state & ACACHE_LOCKBIT)
		*locked = true;
	else
		*locked = false;

	*ncached = state & ACACHE_NCACHED_MASK;
	*low_water = (state >> ACACHE_NCACHED_BITS) & ACACHE_NCACHED_MASK;
	assert(*low_water <= *ncached);

	return state;
}

JEMALLOC_ALWAYS_INLINE bool
cbin_state_commit(arena_cache_bin_t *cbin, const acache_state_t old_state,
    const acache_state_t new_state)
{
	UNUSED unsigned ncached, low_water;

	ncached = new_state & ACACHE_NCACHED_MASK;
	low_water = (new_state >> ACACHE_NCACHED_BITS) & ACACHE_NCACHED_MASK;
	assert(low_water <= ncached);
	assert((uint64_t)(&cbin->data) % sizeof(void *) == 0);

	return atomic_cas_uint64(&cbin->data, old_state, new_state);
}

JEMALLOC_ALWAYS_INLINE acache_state_t
cbin_state_adjust(const acache_state_t state, const int ncached_diff,
    const bool epoch_inc)
{

	if (ncached_diff < 0) {
		assert((int)(state & ACACHE_NCACHED_MASK) >= -ncached_diff);
	}

	return state + ncached_diff + (epoch_inc ? ACACHE_EPOCH_INC : 0);
}

JEMALLOC_ALWAYS_INLINE acache_state_t
cbin_get_lock_state(const acache_state_t state)
{

	assert(!(state & ACACHE_LOCKBIT));
	return cbin_state_adjust(state, 0, true) | ACACHE_LOCKBIT;
}

JEMALLOC_ALWAYS_INLINE bool
cbin_lock(arena_cache_bin_t *cbin, acache_state_t *state)
{
	acache_state_t new_state, old_state;
	bool cas_fail;

	old_state = *state;
	assert(!(old_state & ACACHE_LOCKBIT));

	new_state = cbin_get_lock_state(old_state);
	cas_fail = cbin_state_commit(cbin, old_state, new_state);
	if (!cas_fail)
		*state = new_state;

	return cas_fail;
}

JEMALLOC_ALWAYS_INLINE bool
cbin_lock_and_get_info(arena_cache_bin_t *cbin, acache_state_t *state,
    size_t *ncached, size_t *low_water)
{
	bool bin_locked, cas_fail;
	unsigned retry = 0;
label_retry:
	*state = cbin_state_get(cbin, ncached, low_water, &bin_locked);
	if (bin_locked) {
		return (true);
	}

	/* Try locking the cache bin */
	cas_fail = cbin_lock(cbin, state);
	if (unlikely(cas_fail)) {
		if (retry++ == ACACHE_CAS_FAIL_RETRY_MAX) {
			return (true);
		}
		SPINWAIT_RETRY;
	}

	return (false);
}

JEMALLOC_ALWAYS_INLINE void
cbin_unlock(arena_cache_bin_t *cbin, const acache_state_t state)
{

	ACCESS_ONCE(cbin->data) = state & ~ACACHE_LOCKBIT;
}

JEMALLOC_ALWAYS_INLINE acache_state_t
cbin_state_pack(const acache_state_t state, const size_t ncached,
    const size_t low_water, const bool locked)
{
	acache_state_t new_state;
	uint64_t epoch;

	/* Update epoch when getting the new state. */
	epoch = (state >> ACACHE_EPOCH_OFF) + 1;
	new_state = (epoch << ACACHE_EPOCH_OFF) | ncached |
		(low_water << ACACHE_NCACHED_BITS);

	if (locked) {
		new_state |= ACACHE_LOCKBIT;
	}

	return new_state;
}

/* Return number of items allocated */
JEMALLOC_ALWAYS_INLINE size_t
cbin_alloc_to_tbin(arena_cache_bin_t *cbin, tcache_bin_t *tbin,
    unsigned n_alloc)
{
	acache_state_t state, new_state;
	size_t ncached, low_water, retry;
	bool cas_fail, locked;

	retry = 0;
label_retry:
	state = cbin_state_get(cbin, &ncached, &low_water, &locked);
	if (locked || (ncached == 0))
		return (0);

	if (n_alloc > ncached)
		n_alloc = ncached;

	assert(tbin->ncached == 0);
	/* tbin->avail[-ncached ... -1] are available items. */
	memcpy(&tbin->avail[-(int)n_alloc], &cbin->avail[ncached - n_alloc],
	    sizeof(void *) * n_alloc);

	/* Update state and low_water if necessary */
	if (low_water <= ncached - n_alloc) {
		new_state = cbin_state_adjust(state, -n_alloc, true);
	} else {
		low_water = ncached = ncached - n_alloc;
		new_state = cbin_state_pack(state, ncached, low_water, false);
	}

	cas_fail = cbin_state_commit(cbin, state, new_state);
	if (unlikely(cas_fail)) {
		assert((ACCESS_ONCE(cbin->data) >> ACACHE_EPOCH_OFF) !=
		    (state >> ACACHE_EPOCH_OFF));
		if (retry++ == ACACHE_CAS_FAIL_RETRY_MAX) {
			return (0);
		}
		SPINWAIT_RETRY;
	}
	assert(tbin->avail[-(int)n_alloc] && tbin->avail[-1]);
	tbin->ncached = n_alloc;

	return (n_alloc);
}

/*
 * Allocate a single item from acache if available. Currently used for large
 * allocations only.
 */
JEMALLOC_ALWAYS_INLINE void *
cbin_alloc(arena_cache_bin_t *cbin)
{
	acache_state_t state, new_state;
	size_t ncached, low_water, retry;
	bool cas_fail, locked;
	void *ret;

	retry = 0;
label_retry:
	state = cbin_state_get(cbin, &ncached, &low_water, &locked);
	if (locked || (ncached == 0))
		return (0);

	ret = cbin->avail[ncached - 1];
	/* Update state and low_water if necessary. */
	if (low_water <= ncached - 1) {
		new_state = cbin_state_adjust(state, -1, true);
	} else {
		low_water = ncached = ncached - 1;
		new_state = cbin_state_pack(state, ncached, low_water, false);
	}
	cas_fail = cbin_state_commit(cbin, state, new_state);
	if (unlikely(cas_fail)) {
		assert((ACCESS_ONCE(cbin->data) >> ACACHE_EPOCH_OFF) !=
		    (state >> ACACHE_EPOCH_OFF));
		if (retry++ == ACACHE_CAS_FAIL_RETRY_MAX) {
			return (0);
		}
		SPINWAIT_RETRY;
	}

	return (ret);
}

/* Allocated a single item. Use cached memory if available. */
JEMALLOC_ALWAYS_INLINE void *
arena_cache_alloc_large(tsdn_t *tsdn, arena_t *arena, szind_t binind,
    size_t usize, bool zero)
{
	arena_cache_bin_t *cbin = &arena->acache->cbins[binind];
	void *ret = NULL;

	assert(binind <= nhbins);
	if (!zero) {
		ret = cbin_alloc(cbin);
	}

	return ret ? ret : arena_malloc_large(tsdn, arena, binind, zero);
}

/*
 * Fill tcache from acache if there is memory cached in acache. Return number of
 * items filled.
 */
JEMALLOC_ALWAYS_INLINE size_t
arena_cache_alloc_small(tsdn_t *tsdn, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, szind_t binind)
{
	arena_cache_bin_t *cbin = &arena->acache->cbins[binind];
	size_t nfill;

	assert(binind < NBINS);
	nfill = tcache_bin_info[binind].ncached_max >> tbin->lg_fill_div;
	assert(nfill);

	if (binind < opt_acache_bypass) {
		size_t n_items_last_thd;
		assert(cbin->n_items_last_thd <=
		    arena_cache_bin_info[binind].ncached_max);

		n_items_last_thd = ACCESS_ONCE(cbin->n_items_last_thd);
		/*
		 * Reusing small items (cacheline granularity) from acache may cause serious
		 * fragmentation and false sharing. Only reuse when the current thread was
		 * the one that deallocated the memory into acache.
		 */
		if (cbin->last_thd != tsdn || n_items_last_thd == 0) {
			return 0;
		}

		/* There are race conditions around this, which could cause some inefficient
		 * cache reuse. However it should be rare and still safe. */
		if (n_items_last_thd < nfill) {
			nfill = n_items_last_thd;
		}
		ACCESS_ONCE(cbin->n_items_last_thd) = n_items_last_thd - nfill;
	}

	return cbin_alloc_to_tbin(cbin, tbin, nfill);
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_merge_stats(arena_t *arena, szind_t binind, uint64_t nrequests)
{
	arena_cache_bin_t *cbin = &arena->acache->cbins[binind];
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
	uint64_t queued_nrequests;

	if (is_large) {
		assert(binind >= NBINS);
		arena->stats.nrequests_large += nrequests;
		arena->stats.lstats[binind - NBINS].nrequests += nrequests;
	} else {
		bin->stats.nflushes++;
		bin->stats.nrequests += nrequests;
	}

	if (!config_acache || !opt_acache)
		return;

	/* Called w/ arena or bin lock. Thus read + atomic_sub is safe. */
	queued_nrequests = ACCESS_ONCE(cbin->queued_nrequests);

	if (queued_nrequests) {
		atomic_sub_uint64(&cbin->queued_nrequests, queued_nrequests);

		if (is_large) {
			assert(binind >= NBINS);
			arena->stats.nrequests_large += queued_nrequests;
			arena->stats.lstats[binind - NBINS].nrequests += queued_nrequests;
		} else {
			uint64_t queued_nflushes = ACCESS_ONCE(cbin->queued_nflushes);
			atomic_sub_uint64(&cbin->queued_nflushes, queued_nflushes);
			bin->stats.nflushes += queued_nflushes;
			bin->stats.nrequests += queued_nrequests;
		}
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_flush_locked(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
		uint64_t nrequests, const bool is_large)
{

	if (config_stats) {
		arena_cache_stats_merge_locked(arena, bin, cbin, binind, nrequests,
		    is_large);
	}

	for (size_t i = 0; i < n_items; i++) {
		void *ptr = items[i];
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
		if (is_large) {
			arena_dalloc_large_junked_locked(tsdn, arena, chunk, ptr);
		} else {
			size_t pageind = ((uintptr_t)ptr -
				  (uintptr_t)chunk) >> LG_PAGE;
			arena_chunk_map_bits_t *bitselm =
				  arena_bitselm_get_mutable(chunk, pageind);
			arena_dalloc_bin_junked_locked(tsdn, arena, chunk, ptr, bitselm);
		}
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_flush_small(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
		uint64_t nrequests)
{
	unsigned i, nflush, ndeferred;
	void *ptr;

	for (nflush = n_items; nflush > 0; nflush = ndeferred) {
		/* Process the arena bin associated with the first object. */
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(items[0]);
		arena_t *bin_arena = extent_node_arena_get(&chunk->node);
		/* Needed as acache is enabled / disabled. */
		arena_bin_t *bin = &bin_arena->bins[binind];

		malloc_mutex_lock(tsdn, &bin->lock);
		if (config_stats && bin_arena == arena) {
			arena_cache_stats_merge_locked(arena, bin, cbin, binind, nrequests,
			    false);
		}

		ndeferred = 0;
		for (i = 0; i < nflush; i++) {
			ptr = items[i];
			assert(ptr != NULL);
			chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
			if (extent_node_arena_get(&chunk->node) == bin_arena) {
				size_t pageind = ((uintptr_t)ptr -
				    (uintptr_t)chunk) >> LG_PAGE;
				arena_chunk_map_bits_t *bitselm =
				    arena_bitselm_get_mutable(chunk, pageind);
				arena_dalloc_bin_junked_locked(tsdn, bin_arena, chunk, ptr, bitselm);
			} else {
				/*
				 * This object was allocated via a different
				 * arena bin than the one that is currently
				 * locked.  Stash the object, so that it can be
				 * handled in a future pass.
				 */
				items[ndeferred++] = ptr;
			}
		}

		malloc_mutex_unlock(tsdn, &bin->lock);
		arena_decay_ticks(tsdn, bin_arena, nflush - ndeferred);
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_flush_large(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
		uint64_t nrequests)
{
	unsigned i, nflush, ndeferred;
	void *ptr;

	for (nflush = n_items; nflush > 0; nflush = ndeferred) {
		/* Process the arena associated with the first object. */
		arena_chunk_t *chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(items[0]);
		arena_t *bin_arena = extent_node_arena_get(&chunk->node);

		malloc_mutex_lock(tsdn, &bin_arena->lock);
		if (config_stats && bin_arena == arena) {
			arena_cache_stats_merge_locked(arena, bin, cbin, binind, nrequests,
			   true);
		}

		ndeferred = 0;
		for (i = 0; i < nflush; i++) {
			ptr = items[i];
			assert(ptr != NULL);
			chunk = (arena_chunk_t *)CHUNK_ADDR2BASE(ptr);
			if (extent_node_arena_get(&chunk->node) == bin_arena) {
				arena_dalloc_large_junked_locked(tsdn,
				    bin_arena, chunk, ptr);
			} else {
				/*
				 * This object was allocated via a different
				 * arena than the one that is currently locked.
				 * Stash the object, so that it can be handled
				 * in a future pass.
				 */
				items[ndeferred++] = ptr;
			}
		}
		malloc_mutex_unlock(tsdn, &bin_arena->lock);

		arena_decay_ticks(tsdn, bin_arena, nflush - ndeferred);
	}
}

JEMALLOC_ALWAYS_INLINE void
arena_cache_flush(tsdn_t *tsdn, arena_t *arena, arena_bin_t *bin,
		arena_cache_bin_t *cbin, void **items, size_t n_items, szind_t binind,
		uint64_t nrequests, const bool is_large)
{

	if (!is_large) {
		assert(binind < NBINS);
		arena_cache_flush_small(tsdn, arena, bin, cbin, items, n_items, binind,
	    nrequests);
	} else {
		assert(binind >= NBINS);
		arena_cache_flush_large(tsdn, arena, bin, cbin, items, n_items, binind,
	    nrequests);
	}
}


JEMALLOC_ALWAYS_INLINE void
arena_cache_dalloc(tsdn_t *tsdn, arena_t *arena, void **items,
    size_t n_items, szind_t binind, uint64_t nrequests, const bool is_large)
{
	arena_bin_t *bin = &arena->bins[binind];
	arena_cache_bin_t *cbin = &arena->acache->cbins[binind];
	acache_state_t state;
	size_t ncached, low_water;

	if (!config_acache || !opt_acache ||
	    cbin_lock_and_get_info(cbin, &state, &ncached, &low_water)) {
		/* Flush back to arena if acache is not available. */
		arena_cache_flush(tsdn, arena, bin, cbin, items, n_items, binind,
		    nrequests, is_large);
		return;
	}

	assert(n_items <= arena_cache_bin_info[binind].ncached_max);
	/* Check if need to flush some items back to arena. */
	if (ncached + n_items <= arena_cache_bin_info[binind].ncached_max) {
		/* Queue the stats update as we are not going to lock the arena. */
		arena_cache_merge_stats(arena, binind, nrequests);
	} else {
		if (binind < opt_acache_bypass) {
			assert(!is_large);
			arena_cache_flush_small(tsdn, arena, bin, cbin, cbin->avail, ncached,
			    binind, nrequests);
			state = cbin_state_pack(state, 0, 0, true);
			ncached = low_water = 0;
			cbin->last_thd = NULL;
		} else {
			size_t nkeep, nflush;

			if (n_items >= arena_cache_bin_info[binind].flush_remain) {
				nkeep = 0;
			} else {
				nkeep = arena_cache_bin_info[binind].flush_remain - n_items;
			}
			/*
			 * After flushing back to arena, and return the n_items to acache, we leave
			 * flush_remain items in acache.
			 */
			nflush = ncached - nkeep;
			assert(nkeep < ncached);

			/*
			 * If we use cbin_alloc to get items from the cbin, it's possible to avoid
			 * locking the cbin while flushing to arena. However it requires additional
			 * memcpy in that case.
			 */
			arena_cache_flush(tsdn, arena, bin, cbin, &cbin->avail[nkeep],
			    nflush, binind, nrequests, is_large);

			/* Update state which will be used for unlocking. */
			ncached = nkeep;
			if (low_water > nkeep + n_items) {
				low_water = nkeep + n_items;
				/* ncached will be adjusted when commiting */
				state = cbin_state_pack(state, ncached, low_water, true);
			} else {
				state = cbin_state_adjust(state, -nflush, false);
			}
			/* Fall through to free to cbin and unlock. */
		}
	}
	assert(ncached + n_items <= arena_cache_bin_info[binind].ncached_max);

	/*
	 * Do the memcpy, then release the bitlock & update ncached in a single
	 * store.
	 */
	memcpy(&cbin->avail[ncached], items, sizeof(void *) * n_items);
	cc_barrier(); /* For x86, compiler barrier is sufficient. */
	assert(cbin->avail[ncached + n_items - 1]);
	assert((state & ACACHE_NCACHED_MASK) <=
	    arena_cache_bin_info[binind].ncached_max);

	if (binind < opt_acache_bypass) {
		if (cbin->last_thd == tsdn) {
			cbin->n_items_last_thd += n_items;
			assert(cbin->n_items_last_thd <=
			    arena_cache_bin_info[binind].ncached_max);
		} else {
			cbin->last_thd = tsdn;
			cbin->n_items_last_thd = n_items;
		}
	}

	/* Release the bitlock, and adjust ncached. */
	cbin_unlock(cbin, cbin_state_adjust(state, n_items, false));
}

#  endif /* JEMALLOC_ARENA_INLINE_B */
#endif

#endif /* JEMALLOC_H_INLINES */
/******************************************************************************/
