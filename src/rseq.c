/*
 * rseq.c
 *
 * Copyright (C) 2016 Mathieu Desnoyers <mathieu.desnoyers@efficios.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; only
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 */

//#define _GNU_SOURCE
#include <errno.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <syscall.h>
#include <assert.h>
#include <signal.h>

#include "jemalloc/internal/membarrier.h"
#include "jemalloc/internal/rseq.h"

#ifdef __NR_membarrier
# define membarrier(...)		syscall(__NR_membarrier, __VA_ARGS__)
#else
# define membarrier(...)		-ENOSYS
#endif

struct rseq_thread_state {
	uint32_t fallback_wait_cnt;
	uint32_t fallback_cnt;
	sigset_t sigmask_saved;
};

__attribute__((weak)) __thread volatile struct rseq __rseq_abi = {
	.u.e.cpu_id = -1,
};

static __thread volatile struct rseq_thread_state rseq_thread_state;

int rseq_has_sys_membarrier;

static int sys_rseq(volatile struct rseq *rseq_abi, int flags)
{
  #define __NR_rseq 329
	return syscall(__NR_rseq, rseq_abi, flags);
}

int rseq_register_current_thread(void)
{
	int rc;

	rc = sys_rseq(&__rseq_abi, 0);
	if (rc) {
		fprintf(stderr, "Error: sys_rseq(...) failed(%d): %s\n",
			errno, strerror(errno));
		return -1;
	}
	assert(rseq_current_cpu() >= 0);
	return 0;
}

int rseq_unregister_current_thread(void)
{
	int rc;

	rc = sys_rseq(NULL, 0);
	if (rc) {
		fprintf(stderr, "Error: sys_rseq(...) failed(%d): %s\n",
			errno, strerror(errno));
		return -1;
	}
	return 0;
}

int rseq_init_lock(struct rseq_lock *rlock)
{
	int ret;

	ret = pthread_mutex_init(&rlock->lock, NULL);
	if (ret) {
		errno = ret;
		return -1;
	}
	rlock->state = RSEQ_LOCK_STATE_RESTART;
	return 0;
}

int rseq_destroy_lock(struct rseq_lock *rlock)
{
	int ret;

	ret = pthread_mutex_destroy(&rlock->lock);
	if (ret) {
		errno = ret;
		return -1;
	}
	return 0;
}

static void signal_off_save(sigset_t *oldset)
{
	sigset_t set;
	int ret;

	sigfillset(&set);
	ret = pthread_sigmask(SIG_BLOCK, &set, oldset);
	if (ret)
		abort();
}

static void signal_restore(sigset_t oldset)
{
	int ret;

	ret = pthread_sigmask(SIG_SETMASK, &oldset, NULL);
	if (ret)
		abort();
}

static void rseq_fallback_lock(struct rseq_lock *rlock)
{
	signal_off_save((sigset_t *)&rseq_thread_state.sigmask_saved);
	pthread_mutex_lock(&rlock->lock);
	rseq_thread_state.fallback_cnt++;
	/*
	 * For concurrent threads arriving before we set LOCK:
	 * reading cpu_id after setting the state to LOCK
	 * ensures they restart.
	 */
	ACCESS_ONCE(rlock->state) = RSEQ_LOCK_STATE_LOCK;
	/*
	 * For concurrent threads arriving after we set LOCK:
	 * those will grab the lock, so we are protected by
	 * mutual exclusion.
	 */
}

void rseq_fallback_wait(struct rseq_lock *rlock)
{
	signal_off_save((sigset_t *)&rseq_thread_state.sigmask_saved);
	pthread_mutex_lock(&rlock->lock);
	rseq_thread_state.fallback_wait_cnt++;
	pthread_mutex_unlock(&rlock->lock);
	signal_restore(rseq_thread_state.sigmask_saved);
}

static void rseq_fallback_unlock(struct rseq_lock *rlock, int cpu_at_start)
{
	/*
	 * Concurrent rseq arriving before we set state back to RESTART
	 * grab the lock. Those arriving after we set state back to
	 * RESTART will perform restartable critical sections. The next
	 * owner of the lock will take take of making sure it prevents
	 * concurrent restartable sequences from completing.  We may be
	 * writing from another CPU, so update the state with a store
	 * release semantic to ensure restartable sections will see our
	 * side effect (writing to *p) before they enter their
	 * restartable critical section.
	 *
	 * In cases where we observe that we are on the right CPU after the
	 * critical section, program order ensures that following restartable
	 * critical sections will see our stores, so we don't have to use
	 * store-release or membarrier.
	 *
	 * Use sys_membarrier when available to remove the memory barrier
	 * implied by smp_load_acquire().
	 */
	barrier();
	if (likely(rseq_current_cpu() == cpu_at_start)) {
		ACCESS_ONCE(rlock->state) = RSEQ_LOCK_STATE_RESTART;
	} else {
		if (!has_fast_acquire_release() && rseq_has_sys_membarrier) {
			if (membarrier(MEMBARRIER_CMD_SHARED, 0))
				abort();
			ACCESS_ONCE(rlock->state) = RSEQ_LOCK_STATE_RESTART;
		} else {
			/*
			 * Store with release semantic to ensure
			 * restartable sections will see our side effect
			 * (writing to *p) before they enter their
			 * restartable critical section. Matches
			 * smp_load_acquire() in rseq_start().
			 */
			smp_store_release(&rlock->state,
				RSEQ_LOCK_STATE_RESTART);
		}
	}
	pthread_mutex_unlock(&rlock->lock);
	signal_restore(rseq_thread_state.sigmask_saved);
}

int rseq_fallback_current_cpu(void)
{
	int cpu;

	cpu = sched_getcpu();
	if (cpu < 0) {
		perror("sched_getcpu()");
		abort();
	}
	return cpu;
}

int rseq_fallback_begin(struct rseq_lock *rlock)
{
	rseq_fallback_lock(rlock);
	return rseq_fallback_current_cpu();
}

void rseq_fallback_end(struct rseq_lock *rlock, int cpu)
{
	rseq_fallback_unlock(rlock, cpu);
}

/* Handle non-initialized rseq for this thread. */
void rseq_fallback_noinit(struct rseq_state *rseq_state)
{
	rseq_state->lock_state = RSEQ_LOCK_STATE_FAIL;
	rseq_state->cpu_id = 0;
}

uint32_t rseq_get_fallback_wait_cnt(void)
{
	return rseq_thread_state.fallback_wait_cnt;
}

uint32_t rseq_get_fallback_cnt(void)
{
	return rseq_thread_state.fallback_cnt;
}

void __attribute__((constructor)) rseq_init(void)
{
	int ret;

	ret = membarrier(MEMBARRIER_CMD_QUERY, 0);
	if (ret >= 0 && (ret & MEMBARRIER_CMD_SHARED))
		rseq_has_sys_membarrier = 1;
}
