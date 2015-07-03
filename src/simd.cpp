#include <layer/simd.h>
#if defined(__LINUX__)
#include <malloc.h>
#endif

NAMESPACE_BEGIN(layer)
NAMESPACE_BEGIN(simd)

/* Assumed L1 cache line size for alignment purposes */
#if !defined(L1_CACHE_LINE_SIZE)
#define L1_CACHE_LINE_SIZE 64
#endif

void * malloc(size_t size) {
#if defined(__WINDOWS__)
    return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
#elif defined(__OSX__)
	/* OSX malloc already returns 16-byte aligned data suitable
	   for AltiVec and SSE computations */
    return ::malloc(size);
#else
    return memalign(L1_CACHE_LINE_SIZE, size);
#endif
}

void free(void *ptr) {
#if defined(__WINDOWS__)
	_aligned_free(ptr);
#else
	::free(ptr);
#endif
}

NAMESPACE_END(simd)
NAMESPACE_END(layer)
