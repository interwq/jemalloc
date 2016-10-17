/*
 * rseq-x86.h
 *
 * (C) Copyright 2016 - Mathieu Desnoyers <mathieu.desnoyers@efficios.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifdef __x86_64__

#define smp_mb()	__asm__ __volatile__ ("mfence" : : : "memory")
#define smp_rmb()	barrier()
#define smp_wmb()	barrier()

#define smp_load_acquire(p)						\
__extension__ ({							\
	__typeof(*p) ____p1 = READ_ONCE(*p);				\
	barrier();							\
	____p1;								\
})

#define smp_acquire__after_ctrl_dep()	smp_rmb()

#define smp_store_release(p, v)						\
do {									\
	barrier();							\
	WRITE_ONCE(*p, v);						\
} while (0)

#define has_fast_acquire_release()	1
#define has_single_copy_load_64()	1

/*
 * The __rseq_table section can be used by debuggers to better handle
 * single-stepping through the restartable critical sections.
 */
#define RSEQ_FINISH_ASM(_target_final, _to_write_final, _start_value, \
		_failure, _spec_store, _spec_input, \
		_final_store, _final_input, _extra_clobber, \
		_setup, _teardown, _scratch) \
do { \
	_scratch \
	__asm__ __volatile__ goto ( \
		".pushsection __rseq_table, \"aw\"\n\t" \
		".balign 32\n\t" \
		"3:\n\t" \
		".quad 1f, 2f, 4f, 0x0\n\t" \
		".popsection\n\t" \
		"1:\n\t" \
		_setup \
		RSEQ_INJECT_ASM(1) \
		"movq $3b, %[rseq_cs]\n\t" \
		RSEQ_INJECT_ASM(2) \
		"cmpl %[start_event_counter], %[current_event_counter]\n\t" \
		"jnz 4f\n\t" \
		RSEQ_INJECT_ASM(3) \
		_spec_store \
		_final_store \
		"2:\n\t" \
		RSEQ_INJECT_ASM(5) \
		"movq $0, %[rseq_cs]\n\t" \
		_teardown \
		".pushsection __rseq_failure, \"a\"\n\t" \
		"4:\n\t" \
		"movq $0, %[rseq_cs]\n\t" \
		_teardown \
		"jmp %l[failure]\n\t" \
		".popsection\n\t" \
		: /* gcc asm goto does not allow outputs */ \
		: [start_event_counter]"r"((_start_value).event_counter), \
		  [current_event_counter]"m"((_start_value).rseqp->u.e.event_counter), \
		  [rseq_cs]"m"((_start_value).rseqp->rseq_cs) \
		  _spec_input \
		  _final_input \
		  RSEQ_INJECT_INPUT \
		: "memory", "cc" \
		  _extra_clobber \
		  RSEQ_INJECT_CLOBBER \
		: _failure \
	); \
} while (0)

#define RSEQ_FINISH_FINAL_STORE_ASM() \
		"movq %[to_write_final], %[target_final]\n\t"

/* x86-64 is TSO */
#define RSEQ_FINISH_FINAL_STORE_RELEASE_ASM() \
		RSEQ_FINISH_FINAL_STORE_ASM()

#define RSEQ_FINISH_FINAL_STORE_INPUT(_target_final, _to_write_final) \
		, [to_write_final]"r"(_to_write_final), \
		[target_final]"m"(*(_target_final))

#define RSEQ_FINISH_SPECULATIVE_STORE_ASM() \
		"movq %[to_write_spec], %[target_spec]\n\t" \
		RSEQ_INJECT_ASM(4)

#define RSEQ_FINISH_SPECULATIVE_STORE_INPUT(_target_spec, _to_write_spec) \
		, [to_write_spec]"r"(_to_write_spec), \
		[target_spec]"m"(*(_target_spec))

/* TODO: implement a faster memcpy. */
#define RSEQ_FINISH_MEMCPY_STORE_ASM() \
		"test %[len_memcpy], %[len_memcpy]\n\t" \
		"jz 333f\n\t" \
		"222:\n\t" \
		"movb (%[to_write_memcpy]), %%al\n\t" \
		"movb %%al, (%[target_memcpy])\n\t" \
		"inc %[to_write_memcpy]\n\t" \
		"inc %[target_memcpy]\n\t" \
		"dec %[len_memcpy]\n\t" \
		"jnz 222b\n\t" \
		"333:\n\t" \
		RSEQ_INJECT_ASM(4)

#define RSEQ_FINISH_MEMCPY_STORE_INPUT(_target_memcpy, _to_write_memcpy, _len_memcpy) \
		, [to_write_memcpy]"r"(_to_write_memcpy), \
		[target_memcpy]"r"(_target_memcpy), \
		[len_memcpy]"r"(_len_memcpy), \
		[rseq_scratch0]"m"(rseq_scratch[0]), \
		[rseq_scratch1]"m"(rseq_scratch[1]), \
		[rseq_scratch2]"m"(rseq_scratch[2])

#define RSEQ_FINISH_MEMCPY_CLOBBER()	\
		, "rax"

#define RSEQ_FINISH_MEMCPY_SCRATCH() \
		uint64_t rseq_scratch[3];

/*
 * We need to save and restore those input registers so they can be
 * modified within the assembly.
 */
#define RSEQ_FINISH_MEMCPY_SETUP() \
		"movq %[to_write_memcpy], %[rseq_scratch0]\n\t" \
		"movq %[target_memcpy], %[rseq_scratch1]\n\t" \
		"movq %[len_memcpy], %[rseq_scratch2]\n\t"

#define RSEQ_FINISH_MEMCPY_TEARDOWN() \
		"movq %[rseq_scratch2], %[len_memcpy]\n\t" \
		"movq %[rseq_scratch1], %[target_memcpy]\n\t" \
		"movq %[rseq_scratch0], %[to_write_memcpy]\n\t"

#elif __i386__

/*
 * Support older 32-bit architectures that do not implement fence
 * instructions.
 */
#define smp_mb()	\
	__asm__ __volatile__ ("lock; addl $0,0(%%esp)" : : : "memory")
#define smp_rmb()	\
	__asm__ __volatile__ ("lock; addl $0,0(%%esp)" : : : "memory")
#define smp_wmb()	\
	__asm__ __volatile__ ("lock; addl $0,0(%%esp)" : : : "memory")

#define smp_load_acquire(p)						\
__extension__ ({							\
	__typeof(*p) ____p1 = READ_ONCE(*p);				\
	smp_mb();							\
	____p1;								\
})

#define smp_acquire__after_ctrl_dep()	smp_rmb()

#define smp_store_release(p, v)						\
do {									\
	smp_mb();							\
	WRITE_ONCE(*p, v);						\
} while (0)

#define has_fast_acquire_release()	0
#define has_single_copy_load_64()	0

/*
 * Use eax as scratch register and take memory operands as input to
 * lessen register pressure. Especially needed when compiling
 * do_rseq_memcpy() in O0.
 */
#define RSEQ_FINISH_ASM(_target_final, _to_write_final, _start_value, \
		_failure, _spec_store, _spec_input, \
		_final_store, _final_input, _extra_clobber, \
		_setup, _teardown, _scratch) \
do { \
	_scratch \
	__asm__ __volatile__ goto ( \
		".pushsection __rseq_table, \"aw\"\n\t" \
		".balign 32\n\t" \
		"3:\n\t" \
		".long 1f, 0x0, 2f, 0x0, 4f, 0x0, 0x0, 0x0\n\t" \
		".popsection\n\t" \
		"1:\n\t" \
		_setup \
		RSEQ_INJECT_ASM(1) \
		"movl $3b, %[rseq_cs]\n\t" \
		RSEQ_INJECT_ASM(2) \
		"movl %[start_event_counter], %%eax\n\t" \
		"cmpl %%eax, %[current_event_counter]\n\t" \
		"jnz 4f\n\t" \
		RSEQ_INJECT_ASM(3) \
		_spec_store \
		_final_store \
		"2:\n\t" \
		RSEQ_INJECT_ASM(5) \
		"movl $0, %[rseq_cs]\n\t" \
		_teardown \
		".pushsection __rseq_failure, \"a\"\n\t" \
		"4:\n\t" \
		"movl $0, %[rseq_cs]\n\t" \
		_teardown \
		"jmp %l[failure]\n\t" \
		".popsection\n\t" \
		: /* gcc asm goto does not allow outputs */ \
		: [start_event_counter]"m"((_start_value).event_counter), \
		  [current_event_counter]"m"((_start_value).rseqp->u.e.event_counter), \
		  [rseq_cs]"m"((_start_value).rseqp->rseq_cs) \
		  _spec_input \
		  _final_input \
		  RSEQ_INJECT_INPUT \
		: "memory", "cc", "eax" \
		  _extra_clobber \
		  RSEQ_INJECT_CLOBBER \
		: _failure \
	); \
} while (0)

#define RSEQ_FINISH_FINAL_STORE_ASM() \
		"movl %[to_write_final], %%eax\n\t" \
		"movl %%eax, %[target_final]\n\t"

#define RSEQ_FINISH_FINAL_STORE_RELEASE_ASM() \
		"lock; addl $0,0(%%esp)\n\t" \
		RSEQ_FINISH_FINAL_STORE_ASM()

#define RSEQ_FINISH_FINAL_STORE_INPUT(_target_final, _to_write_final) \
		, [to_write_final]"m"(_to_write_final), \
		[target_final]"m"(*(_target_final))

#define RSEQ_FINISH_SPECULATIVE_STORE_ASM() \
		"movl %[to_write_spec], %%eax\n\t" \
		"movl %%eax, %[target_spec]\n\t" \
		RSEQ_INJECT_ASM(4)

#define RSEQ_FINISH_SPECULATIVE_STORE_INPUT(_target_spec, _to_write_spec) \
		, [to_write_spec]"m"(_to_write_spec), \
		[target_spec]"m"(*(_target_spec))

/* TODO: implement a faster memcpy. */
#define RSEQ_FINISH_MEMCPY_STORE_ASM() \
		"movl %[len_memcpy], %%eax\n\t" \
		"test %%eax, %%eax\n\t" \
		"jz 333f\n\t" \
		"222:\n\t" \
		"movb (%[to_write_memcpy]), %%al\n\t" \
		"movb %%al, (%[target_memcpy])\n\t" \
		"inc %[to_write_memcpy]\n\t" \
		"inc %[target_memcpy]\n\t" \
		"decl %[rseq_scratch2]\n\t" \
		"jnz 222b\n\t" \
		"333:\n\t" \
		RSEQ_INJECT_ASM(4)

#define RSEQ_FINISH_MEMCPY_STORE_INPUT(_target_memcpy, _to_write_memcpy, _len_memcpy) \
		, [to_write_memcpy]"r"(_to_write_memcpy), \
		[target_memcpy]"r"(_target_memcpy), \
		[len_memcpy]"m"(_len_memcpy), \
		[rseq_scratch0]"m"(rseq_scratch[0]), \
		[rseq_scratch1]"m"(rseq_scratch[1]), \
		[rseq_scratch2]"m"(rseq_scratch[2])

#define RSEQ_FINISH_MEMCPY_CLOBBER()

#define RSEQ_FINISH_MEMCPY_SCRATCH() \
		uint32_t rseq_scratch[3];

/*
 * We need to save and restore those input registers so they can be
 * modified within the assembly.
 */
#define RSEQ_FINISH_MEMCPY_SETUP() \
		"movl %[to_write_memcpy], %[rseq_scratch0]\n\t" \
		"movl %[target_memcpy], %[rseq_scratch1]\n\t" \
		"movl %[len_memcpy], %%eax\n\t" \
		"movl %%eax, %[rseq_scratch2]\n\t"

#define RSEQ_FINISH_MEMCPY_TEARDOWN() \
		"movl %[rseq_scratch1], %[target_memcpy]\n\t" \
		"movl %[rseq_scratch0], %[to_write_memcpy]\n\t"

#endif
