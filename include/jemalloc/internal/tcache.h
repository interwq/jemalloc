/******************************************************************************/
#ifdef JEMALLOC_H_TYPES

typedef struct ccache_bin_s ccache_bin_t;
typedef struct ccache_s ccache_t;
typedef struct tcache_bin_info_s tcache_bin_info_t;
typedef struct tcache_bin_s tcache_bin_t;
typedef struct tcache_s tcache_t;
typedef struct tcaches_s tcaches_t;

/*
 * tcache pointers close to NULL are used to encode state information that is
 * used for two purposes: preventing thread caching on a per thread basis and
 * cleaning up during thread shutdown.
 */
#define	TCACHE_STATE_DISABLED		((tcache_t *)(uintptr_t)1)
#define	TCACHE_STATE_REINCARNATED	((tcache_t *)(uintptr_t)2)
#define	TCACHE_STATE_PURGATORY		((tcache_t *)(uintptr_t)3)
#define	TCACHE_STATE_MAX		TCACHE_STATE_PURGATORY

/*
 * Absolute minimum number of cache slots for each small bin.
 */
#define	TCACHE_NSLOTS_SMALL_MIN		20

/*
 * Absolute maximum number of cache slots for each small bin in the thread
 * cache.  This is an additional constraint beyond that imposed as: twice the
 * number of regions per run for this size class.
 *
 * This constant must be an even number.
 */
#define	TCACHE_NSLOTS_SMALL_MAX		200

/* Number of cache slots for large size classes. */
#define	TCACHE_NSLOTS_LARGE		20

/* (1U << opt_lg_tcache_max) is used to compute tcache_maxclass. */
#define	LG_TCACHE_MAXCLASS_DEFAULT	15

/*
 * TCACHE_GC_SWEEP is the approximate number of allocation events between
 * full GC sweeps.  Integer rounding may cause the actual number to be
 * slightly higher, since GC is performed incrementally.
 */
#define	TCACHE_GC_SWEEP			(8192)

/* Number of tcache allocation/deallocation events between incremental GCs. */
#define	TCACHE_GC_INCR							\
    ((TCACHE_GC_SWEEP / NBINS) + ((TCACHE_GC_SWEEP / NBINS == 0) ? 0 : 1))

#endif /* JEMALLOC_H_TYPES */
/******************************************************************************/
#ifdef JEMALLOC_H_STRUCTS

typedef enum {
	tcache_enabled_false   = 0, /* Enable cast to/from bool. */
	tcache_enabled_true    = 1,
	tcache_enabled_default = 2
} tcache_enabled_t;

/*
 * Read-only information associated with each element of tcache_t's tbins array
 * is stored separately, mainly to reduce memory usage.
 */
struct tcache_bin_info_s {
	unsigned	ncached_max;	/* Upper limit on ncached. */
};

/* struct tcache_bin_s { */
/* 	tcache_bin_stats_t tstats; */
/* 	unsigned	lg_fill_div;	/\* Fill (ncached_max >> lg_fill_div). *\/ */

/* 	int		low_water;	/\* Min # cached since last GC. *\/ */
/* 	unsigned	ncached;	/\* # of cached objects. *\/ */

/* 	/\* */
/* 	 * To make use of adjacent cacheline prefetch, the items in the avail */
/* 	 * stack goes to higher address for newer allocations.  avail points */
/* 	 * just above the available space, which means that */
/* 	 * avail[-ncached, ... -1] are available items and the lowest item will */
/* 	 * be allocated first. */
/* 	 *\/ */
/* 	void		**avail;	/\* Stack of available objects. *\/ */
/* }; */

typedef enum {
	ccache_success	= 0, /* Common case */
    ccache_empty	= 1,
	ccache_retry	= 2
} ret_status_t;

struct ccache_bin_s {
    //FIXME: deal with 32bit cpus? need this be a single word.
    /*
     *  last thread (32bits) + lock (1bit) + low_water (10bits) + ncached lowest 10bits
     */
    /* The lock bit can only be modified with cpu_lock. */
    uint64_t    data;
//	int		low_water;	/* Min # cached since last GC. */
//	unsigned	ncached;	/* # of cached objects. */
    bool refilled;

	unsigned	lg_fill_div;	/* Fill (ncached_max >> lg_fill_div). */

	/*
	 * To make use of adjacent cacheline prefetch, the items in the avail
	 * stack goes to higher address for newer allocations.  avail points
	 * just above the available space, which means that
	 * avail[-ncached, ... -1] are available items and the lowest item will
	 * be allocated first.
	 */
	void		**avail;	/* Stack of available objects. */

    // this won't be touched on the fastest path (cached allocation)
    malloc_mutex_t cpu_lock;
};

struct ccache_s {
    unsigned		ind;
    szind_t		next_gc_bin;	/* Next bin to GC. */
	ccache_bin_t	cbins[1];	/* Dynamically sized. */
};

struct tcache_bin_s {
    // if not using percpu caching, add the actual bins here.
	tcache_bin_stats_t tstats;
};

struct tcache_s {
	ql_elm(tcache_t) link;		 /* Used for aggregating stats. */
	uint64_t	prof_accumbytes; /* Cleared after arena_prof_accum(). */
	ticker_t	gc_ticker;       /* Drives incremental GC. */
//	szind_t		next_gc_bin;	/* Next bin to GC. */

	tcache_bin_t	tbins[1];
	/*
	 * The pointer stacks associated with tbins follow as a contiguous
	 * array.  During tcache initialization, the avail pointer in each
	 * element of tbins is initialized to point to the proper offset within
	 * this array.
	 */
};

/* Linkage for list of available (previously used) explicit tcache IDs. */
struct tcaches_s {
	union {
		tcache_t	*tcache;
		tcaches_t	*next;
	};
};

#endif /* JEMALLOC_H_STRUCTS */
/******************************************************************************/
#ifdef JEMALLOC_H_EXTERNS

extern bool	opt_tcache;
extern ssize_t	opt_lg_tcache_max;

extern tcache_bin_info_t	*tcache_bin_info;

/*
 * Number of tcache bins.  There are NBINS small-object bins, plus 0 or more
 * large-object bins.
 */
extern unsigned	nhbins;

/* Maximum cached size class. */
extern size_t	tcache_maxclass;

/*
 * Explicit tcaches, managed via the tcache.{create,flush,destroy} mallctls and
 * usable via the MALLOCX_TCACHE() flag.  The automatic per thread tcaches are
 * completely disjoint from this data structure.  tcaches starts off as a sparse
 * array, so it has no physical memory footprint until individual pages are
 * touched.  This allows the entire array to be allocated the first time an
 * explicit tcache is created without a disproportionate impact on memory usage.
 */
extern tcaches_t	*tcaches;

/* TODO */
extern ccache_t	**ccaches;

JEMALLOC_ALWAYS_INLINE unsigned long malloc_getcpu(void);
size_t	tcache_salloc(const void *ptr);
void	tcache_event_hard(tsd_t *tsd, tcache_t *tcache);
void	*tcache_alloc_small_hard(tsd_t *tsd, arena_t *arena, tcache_t *tcache,
    tcache_bin_t *tbin, ccache_bin_t *cbin, szind_t binind, uint32_t curr_thd,
    ret_status_t *ret_status);
uint64_t tcache_bin_flush_small_locked(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem);
void	tcache_bin_flush_small(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem);
uint64_t tcache_bin_flush_large_locked(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind, unsigned rem);
void	tcache_bin_flush_large(tsd_t *tsd, tcache_t *tcache, tcache_bin_t *tbin,
    ccache_bin_t *cbin, szind_t binind,unsigned rem);
void	tcache_arena_associate(tcache_t *tcache, arena_t *arena);
void	tcache_arena_reassociate(tcache_t *tcache, arena_t *oldarena,
    arena_t *newarena);
void	tcache_arena_dissociate(tcache_t *tcache, arena_t *arena);
tcache_t *tcache_get_hard(tsd_t *tsd);
tcache_t *tcache_create(tsd_t *tsd, arena_t *arena);
void	tcache_cleanup(tsd_t *tsd);
void	tcache_enabled_cleanup(tsd_t *tsd);
void	tcache_stats_merge(tcache_t *tcache, arena_t *arena);
bool	ccaches_create(unsigned ncaches);
bool	tcaches_create(tsd_t *tsd, unsigned *r_ind);
void	tcaches_flush(tsd_t *tsd, unsigned ind);
void	tcaches_destroy(tsd_t *tsd, unsigned ind);
bool	tcache_boot(void);

#endif /* JEMALLOC_H_EXTERNS */
/******************************************************************************/
#ifdef JEMALLOC_H_INLINES

#ifndef JEMALLOC_ENABLE_INLINE
void	tcache_event(tsd_t *tsd, tcache_t *tcache);
void	tcache_flush(void);
bool	tcache_enabled_get(void);
tcache_t *tcache_get(tsd_t *tsd, bool create);
void	tcache_enabled_set(bool enabled);
void	*tcache_alloc_easy(ccache_bin_t *cbin, uint32_t curr_thd, ret_status_t *ret_status);//(tcache_bin_t *tbin, bool *ret_status);
void	*tcache_alloc_small(tsd_t *tsd, arena_t *arena, tcache_t *tcache,
    size_t size, szind_t ind, bool zero, bool slow_path);
void	*tcache_alloc_large(tsd_t *tsd, arena_t *arena, tcache_t *tcache,
    size_t size, szind_t ind, bool zero, bool slow_path);
void	tcache_dalloc_small(tsd_t *tsd, tcache_t *tcache, void *ptr,
    szind_t binind, bool slow_path);
void	tcache_dalloc_large(tsd_t *tsd, tcache_t *tcache, void *ptr,
    size_t size, bool slow_path);
tcache_t	*tcaches_get(tsd_t *tsd, unsigned ind);
#endif

#if (defined(JEMALLOC_ENABLE_INLINE) || defined(JEMALLOC_TCACHE_C_))
JEMALLOC_INLINE void
tcache_flush(void)
{
	tsd_t *tsd;

	cassert(config_tcache);

	tsd = tsd_fetch();
	tcache_cleanup(tsd);
}

JEMALLOC_INLINE bool
tcache_enabled_get(void)
{
	tsd_t *tsd;
	tcache_enabled_t tcache_enabled;

	cassert(config_tcache);

	tsd = tsd_fetch();
	tcache_enabled = tsd_tcache_enabled_get(tsd);
	if (tcache_enabled == tcache_enabled_default) {
		tcache_enabled = (tcache_enabled_t)opt_tcache;
		tsd_tcache_enabled_set(tsd, tcache_enabled);
	}

	return ((bool)tcache_enabled);
}

JEMALLOC_INLINE void
tcache_enabled_set(bool enabled)
{
	tsd_t *tsd;
	tcache_enabled_t tcache_enabled;

	cassert(config_tcache);

	tsd = tsd_fetch();

	tcache_enabled = (tcache_enabled_t)enabled;
	tsd_tcache_enabled_set(tsd, tcache_enabled);

	if (!enabled)
		tcache_cleanup(tsd);
}

JEMALLOC_ALWAYS_INLINE tcache_t *
tcache_get(tsd_t *tsd, bool create)
{
	tcache_t *tcache;
    unsigned cpu;

	if (!config_tcache)
		return (NULL);

	tcache = tsd_tcache_get(tsd);
	if (!create)
		return (tcache);
	if (unlikely(tcache == NULL) && tsd_nominal(tsd)) {
		tcache = tcache_get_hard(tsd);
		tsd_tcache_set(tsd, tcache);
        /* Initialize info needed for PerCPU cache. */
        cpu = malloc_getcpu();
        tsd_ccache_set(tsd, ccaches[cpu]);
        // TODO: compatability check, make sure uint32 is good for all platforms.
        tsd_tid_set(tsd, (uint32_t)syscall(SYS_gettid));
	}

	return (tcache);
}

JEMALLOC_ALWAYS_INLINE void
tcache_event(tsd_t *tsd, tcache_t *tcache)
{

	if (TCACHE_GC_INCR == 0)
		return;

	if (unlikely(ticker_tick(&tcache->gc_ticker)))
		tcache_event_hard(tsd, tcache);
}

#define ACCESS_ONCE(x) (*(volatile typeof(x) *)&(x))

#define NCACHED_BITS  10
#define LOW_WATER_BITS NCACHED_BITS
#define TID_BITSOFF 32
#define CBIN_DATA_MASK ((1ULL << TID_BITSOFF) - 1)

#define CBIN_BITLOCK (1UL << 31)
#define CCACHE_IS_LOCKED(v) ((v & CBIN_BITLOCK))

JEMALLOC_ALWAYS_INLINE void
cache_bin_get_info(uint64_t val, uint32_t *ncached, uint32_t *low_water,
    uint32_t *last_thd)
{
    // TODO: pass cbin instead of val into this.
    // last_thd helps us detect thread migration, and avoid ABA problem.
    *last_thd = val >> TID_BITSOFF;
    *ncached = val & ((1UL << NCACHED_BITS) - 1);
    val >>= NCACHED_BITS;
    *low_water = val & ((1UL << LOW_WATER_BITS) - 1);
}

JEMALLOC_ALWAYS_INLINE uint64_t
cbin_data_pack(uint32_t ncached, uint32_t low_water,uint32_t tid)
{
    return ((uint64_t)tid << TID_BITSOFF) | (low_water << NCACHED_BITS)
        | ncached;
}

/* Take the percpu mutex, and set the CBIN_BITLOCK */
JEMALLOC_ALWAYS_INLINE void
ccache_bin_lock(tsd_t *tsd, ccache_bin_t *cbin)
{
    uint64_t new_val, val;
    bool cas_fail;

    malloc_mutex_lock(&cbin->cpu_lock);
retry:
    /* May still be locked by free() */
    while (CCACHE_IS_LOCKED((val = ACCESS_ONCE(cbin->data)))) {
        /* TODO: */
        // Though no guarantee, still be better to use
        // pthread_cond_timedwait();
        struct timespec t;

//        printf("locked bin on cpu %u\n", tsd_ccache_get(tsd)->ind);
        t.tv_sec = 0;
        t.tv_nsec = 100000;
        nanosleep(&t, NULL);
    }

    new_val = ((uint64_t)tsd_tid_get(tsd) << TID_BITSOFF) |
        ((val & CBIN_DATA_MASK) | CBIN_BITLOCK);
    cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
    if (unlikely(cas_fail)) {
        /* The cas can fail (alloc / free may still proceed before we set the
         * lock bit), but the bin won't be locked by other threads. */
        goto retry;
    }
}

/* Clear the CBIN_BITLOCK, then release the percpu mutex. */
JEMALLOC_ALWAYS_INLINE void
ccache_bin_unlock(tsd_t *tsd, ccache_bin_t *cbin, uint64_t new_val)
{
    UNUSED bool cas_fail;
    uint64_t val = ACCESS_ONCE(cbin->data);

    assert(CCACHE_IS_LOCKED(val));
    /* We must be the current owner of the bin. */
    assert(val >> TID_BITSOFF == tsd_tid_get(tsd));
    assert(!(new_val & CBIN_BITLOCK));
    cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
    assert(!cas_fail);
    malloc_mutex_unlock(&cbin->cpu_lock);
}

JEMALLOC_ALWAYS_INLINE void *
tcache_alloc_easy(ccache_bin_t *cbin, uint32_t curr_thd, ret_status_t *ret_status)
{
	void *ret;
    bool cas_fail;
    uint32_t ncached, low_water, last_thd;
    uint64_t val, new_val;

    val = ACCESS_ONCE(cbin->data);
    cache_bin_get_info(val, &ncached, &low_water, &last_thd);
    if (unlikely(CCACHE_IS_LOCKED(val) || last_thd != curr_thd)) {
        *ret_status = ccache_retry;
        return (NULL);
    }
	if (unlikely(ncached == 0)) {
        cbin->refilled = true;
		*ret_status = ccache_empty;
		return (NULL);
	}

	ret = *(cbin->avail - ncached);
    if (unlikely(ncached < low_water)) {
        new_val = cbin_data_pack(ncached - 1, ncached - 1, curr_thd);
    } else {
        // ncached is in lowest 10 bits.
        new_val = val - 1;
    }

    /*
     * Reading ret needs to be done before the atomic cas. No memory barrier
     * needed for x86 (because of the cas) -- for other arch, possibly we need
     * some barrier here).
     */
    cas_fail = atomic_cas_uint64(&cbin->data, val, new_val);
    if (cas_fail) {
        /* Race happened. Some other thread must have modified the owner tid of
         * the bin. */
        assert((ACCESS_ONCE(cbin->data) >> TID_BITSOFF) != curr_thd);
        *ret_status = ccache_retry;
		return (NULL);
    }
	/*
	 * ret_status (instead of ret) should be checked upon the return of
	 * this function.  We avoid checking (ret == NULL) because there is
	 * never a null stored on the avail stack (which is unknown to the
	 * compiler), and eagerly checking ret would cause pipeline stall
	 * (waiting for the cacheline).
	 */
	*ret_status = ccache_success;

	return (ret);
}

JEMALLOC_ALWAYS_INLINE unsigned long long native_read_tscp(unsigned int *aux)
{
    unsigned long low, high;
    asm volatile(".byte 0x0f,0x01,0xf9": "=a" (low), "=d" (high), "=c" (*aux));
    return low | ((uint64_t)high << 32);
}

JEMALLOC_ALWAYS_INLINE unsigned long tscp_getcpu(void)
{
    unsigned long low, high, aux, ret ;
    asm volatile(".byte 0x0f,0x01,0xf9": "=a" (low), "=d" (high), "=c" (aux));
    ret = aux & ((1U << 12) - 1);
//    assert(ret == sched_getcpu());
    return ret;
}

JEMALLOC_ALWAYS_INLINE unsigned long malloc_getcpu(void)
{
    return tscp_getcpu();
}

JEMALLOC_ALWAYS_INLINE void percpu_cache_arena_update(tsd_t *tsd, uint32_t cpu)
{
    ccache_t *ccache = ccaches[cpu];

    if (unlikely(tsd_ccache_get(tsd) != ccache)) {
        /* This should only happen after a thread migration. */
        tsd_ccache_set(tsd, ccache);
    }
#ifdef PERCPU_ARENA
    arena_t *oldarena;
    unsigned oldind;

    oldarena = tsd_arena_get(tsd);
    assert(oldarena);
    oldind = oldarena->ind;
    assert(oldarena == arenas[oldind]);

    if (unlikely(oldind != cpu)) {
        unsigned newind = cpu;
        arena_t *newarena = arena_get(newind, true);
        assert(newarena);

		/* Set new arena/tcache associations. */
		arena_migrate(tsd, oldind, newind);
		if (config_tcache) {
			tcache_t *tcache = tsd_tcache_get(tsd);
			assert (tcache != NULL);
            tcache_arena_reassociate(tcache, oldarena,
				newarena);
		}
    }
#endif
}

/* Read the current CPU id and update current thread's ccache / arena
 * association. Also marks the current thread in the required cbin. */
JEMALLOC_ALWAYS_INLINE void
percpu_info_sync(tsd_t *tsd, szind_t binind)
{
    ccache_t *ccache;
    ccache_bin_t *cbin;

    uint64_t new_val, val;
    uint32_t cpu, prev_thd, curr_thd;
    bool cas_fail;
retry:
    cpu = malloc_getcpu();
    ccache = ccaches[cpu];
    cbin = &ccache->cbins[binind];
    percpu_cache_arena_update(tsd, cpu);

    val = ACCESS_ONCE(cbin->data);
    if (CCACHE_IS_LOCKED(val)) {
        {
//            if (cbin->data
//            printf("observe locked bin on cpu %u\n", cpu);
        }
        /* If the bin is locked, wait for it. */
        malloc_mutex_lock(&cbin->cpu_lock);
        val = ACCESS_ONCE(cbin->data);
        if (CCACHE_IS_LOCKED(val)) {
            struct timespec t;
//            printf("locked bin on cpu %u\n", tsd_ccache_get(tsd)->ind);

            t.tv_sec = 0;
            t.tv_nsec = 100000;
            nanosleep(&t, NULL);

            // pthread_cond_timedwait();
        }

        malloc_mutex_unlock(&cbin->cpu_lock);
        goto retry;
    }

    curr_thd = tsd_tid_get(tsd);
    prev_thd = val >> TID_BITSOFF;
    if (curr_thd == prev_thd)
        return;

    new_val = ((uint64_t)curr_thd << TID_BITSOFF) | (val & CBIN_DATA_MASK);
    if ((cas_fail = atomic_cas_uint64(&cbin->data, val, new_val)))
        goto retry;
}

JEMALLOC_ALWAYS_INLINE void *
tcache_alloc_small(tsd_t *tsd, arena_t *arena, tcache_t *tcache, size_t size,
    szind_t binind, bool zero, bool slow_path)
{
	void *ret;
    ccache_bin_t *cbin;
    tcache_bin_t *tbin;
	ret_status_t ret_status;
	size_t usize JEMALLOC_CC_SILENCE_INIT(0);
    uint32_t curr_thd;

	assert(binind < NBINS);
    curr_thd = tsd_tid_get(tsd);
    tbin = &tcache->tbins[binind];
retry:
    cbin = &tsd_ccache_get(tsd)->cbins[binind];
	ret = tcache_alloc_easy(cbin, curr_thd, &ret_status);
    assert((ret_status == ccache_success) == (ret != NULL));
	if (unlikely(ret_status != ccache_success)) {
        if (ret_status == ccache_retry) {
            /* Thread migration or switch happened. */
            percpu_info_sync(tsd, binind);
            goto retry;
        }

        assert(ret_status == ccache_empty);
		arena = arena_choose(tsd, arena);
		if (unlikely(arena == NULL))
			return (NULL);

		ret = tcache_alloc_small_hard(tsd, arena, tcache, tbin, cbin, binind,
            curr_thd, &ret_status);
        if (ret_status == ccache_retry)
            goto retry;
		if (ret_status != ccache_success) {
            assert(ret == NULL && ret_status == ccache_empty);
            /* If no preemption happened and nothing was returned, we are out of
             * memory. */
			return (NULL);
        }
	}

	assert(ret);
	/*
	 * Only compute usize if required.  The checks in the following if
	 * statement are all static.
	 */
	if (config_prof || (slow_path && config_fill) || unlikely(zero)) {
		usize = index2size(binind);
		assert(tcache_salloc(ret) == usize);
	}

	if (likely(!zero)) {
		if (slow_path && config_fill) {
			if (unlikely(opt_junk_alloc)) {
				arena_alloc_junk_small(ret,
				    &arena_bin_info[binind], false);
			} else if (unlikely(opt_zero))
				memset(ret, 0, usize);
		}
	} else {
		if (slow_path && config_fill && unlikely(opt_junk_alloc)) {
			arena_alloc_junk_small(ret, &arena_bin_info[binind], true);
		}
		memset(ret, 0, usize);
	}

	if (config_stats)
		tbin->tstats.nrequests++;
	if (config_prof)
		tcache->prof_accumbytes += usize;
	tcache_event(tsd, tcache);
	return (ret);
}

JEMALLOC_ALWAYS_INLINE void *
tcache_alloc_large(tsd_t *tsd, arena_t *arena, tcache_t *tcache, size_t size,
    szind_t binind, bool zero, bool slow_path)
{
    ccache_bin_t *cbin;
	void *ret;
	ret_status_t ret_status;
    uint32_t curr_thd;

	assert(binind < nhbins);
    curr_thd = tsd_tid_get(tsd);
retry:
    cbin = &tsd_ccache_get(tsd)->cbins[binind];
	ret = tcache_alloc_easy(cbin, curr_thd, &ret_status);
	assert((ret_status == ccache_success) == (ret != NULL));
    if (unlikely(ret_status != ccache_success)) {
        if (ret_status == ccache_retry) {
            /* Thread migration or switch happened. */
            percpu_info_sync(tsd, binind);
            goto retry;
        }
		/*
		 * Only allocate one large object at a time, because it's quite
		 * expensive to create one and not use it.
		 */
		arena = arena_choose(tsd, arena);
		if (unlikely(arena == NULL))
			return (NULL);

		ret = arena_malloc_large(tsd, arena, binind, zero);
		if (ret == NULL)
			return (NULL);
	} else {
		size_t usize JEMALLOC_CC_SILENCE_INIT(0);

		/* Only compute usize on demand */
		if (config_prof || (slow_path && config_fill) ||
		    unlikely(zero)) {
			usize = index2size(binind);
			assert(usize <= tcache_maxclass);
		}

		if (config_prof && usize == LARGE_MINCLASS) {
			arena_chunk_t *chunk =
			    (arena_chunk_t *)CHUNK_ADDR2BASE(ret);
			size_t pageind = (((uintptr_t)ret - (uintptr_t)chunk) >>
			    LG_PAGE);
			arena_mapbits_large_binind_set(chunk, pageind,
			    BININD_INVALID);
		}
		if (likely(!zero)) {
			if (slow_path && config_fill) {
				if (unlikely(opt_junk_alloc))
					memset(ret, 0xa5, usize);
				else if (unlikely(opt_zero))
					memset(ret, 0, usize);
			}
		} else
			memset(ret, 0, usize);

		if (config_stats)
			tcache->tbins[binind].tstats.nrequests++;
		if (config_prof)
			tcache->prof_accumbytes += usize;
	}

	tcache_event(tsd, tcache);
	return (ret);
}

/* Compiler barrier. */
static inline void
barrier(void)
{
    asm volatile ("" ::: "memory");
}

JEMALLOC_ALWAYS_INLINE void
tcache_dalloc_small(tsd_t *tsd, tcache_t *tcache, void *ptr, szind_t binind,
    bool slow_path)
{
    bool cas_fail;
    ccache_bin_t *cbin;
	tcache_bin_t *tbin;
	tcache_bin_info_t *tbin_info;
    uint32_t ncached, last_thd, curr_thd, low_water;
    uint64_t val;

	assert(tcache_salloc(ptr) <= SMALL_MAXCLASS);
	if (slow_path && config_fill && unlikely(opt_junk_free))
		arena_dalloc_junk_small(ptr, &arena_bin_info[binind]);

    curr_thd = tsd_tid_get(tsd);
    tbin = &tcache->tbins[binind];
	tbin_info = &tcache_bin_info[binind];
retry:
    cbin = &tsd_ccache_get(tsd)->cbins[binind];
    val = ACCESS_ONCE(cbin->data);
    cache_bin_get_info(val, &ncached, &low_water, &last_thd);
    if (unlikely(CCACHE_IS_LOCKED(val) || last_thd != curr_thd)) {
        percpu_info_sync(tsd, binind);
        goto retry;
    }
	if (unlikely(ncached == tbin_info->ncached_max)) {
		tcache_bin_flush_small(tsd, tcache, tbin, cbin, binind,
		    (tbin_info->ncached_max >> 1));
        goto retry;
	}

    /* Lock the bin -- so that we can do the 2 stores next. */
    cas_fail = atomic_cas_uint64(&cbin->data, val, val | CBIN_BITLOCK);
    if (cas_fail)
        goto retry;

    /* We should release the bitlock ASAP. So keep it short. */
	*(cbin->avail - ncached - 1) = ptr;
    barrier(); /* For x86, compiler barrier is sufficient. */
    /* Release the bitlock, and increment ncached (in lowest 10 bits). */
    cbin->data = val + 1;

    barrier(); /* Compiler barrier to make sure shortest critical section. */

    /* Make the stores visible before checking. */
    /* mb_write(); // EXPENSIVE! */
    /* if (unlikely(cbin->has_lock_waiting)) { */
    /*     signal()  */
    /* } */

	tcache_event(tsd, tcache);
}

JEMALLOC_ALWAYS_INLINE void
tcache_dalloc_large(tsd_t *tsd, tcache_t *tcache, void *ptr, size_t size,
    bool slow_path)
{
    bool cas_fail;
	szind_t binind;
    ccache_bin_t *cbin;
	tcache_bin_t *tbin;
	tcache_bin_info_t *tbin_info;
    uint32_t ncached, last_thd, curr_thd, low_water;
    uint64_t val;

	assert((size & PAGE_MASK) == 0);
	assert(tcache_salloc(ptr) > SMALL_MAXCLASS);
	assert(tcache_salloc(ptr) <= tcache_maxclass);

    curr_thd = tsd_tid_get(tsd);
	binind = size2index(size);
	if (slow_path && config_fill && unlikely(opt_junk_free))
		arena_dalloc_junk_large(ptr, size);

	tbin = &tcache->tbins[binind];
	tbin_info = &tcache_bin_info[binind];
retry:
    cbin = &tsd_ccache_get(tsd)->cbins[binind];
    val = ACCESS_ONCE(cbin->data);
    cache_bin_get_info(val, &ncached, &low_water, &last_thd);
    if (unlikely(CCACHE_IS_LOCKED(val) || last_thd != curr_thd)) {
        percpu_info_sync(tsd, binind);
        goto retry;
    }
    assert(ncached <= tbin_info->ncached_max);
	if (unlikely(ncached == tbin_info->ncached_max)) {
		tcache_bin_flush_large(tsd, tcache, tbin, cbin, binind,
		    (tbin_info->ncached_max >> 1));
        goto retry;
	}

    cas_fail = atomic_cas_uint64(&cbin->data, val, val | CBIN_BITLOCK);
    if (cas_fail)
        goto retry;

    /* We should release the bitlock ASAP. So keep it short. */
	*(cbin->avail - ncached - 1) = ptr;
    barrier(); /* For x86, compiler barrier is sufficient. */
    /* Release the bitlock, and increment ncached (in lowest 10 bits). */
    cbin->data = val + 1;

	tcache_event(tsd, tcache);
}

JEMALLOC_ALWAYS_INLINE tcache_t *
tcaches_get(tsd_t *tsd, unsigned ind)
{
	tcaches_t *elm = &tcaches[ind];
	if (unlikely(elm->tcache == NULL))
		elm->tcache = tcache_create(tsd, arena_choose(tsd, NULL));
	return (elm->tcache);
}
#endif

#endif /* JEMALLOC_H_INLINES */
/******************************************************************************/
