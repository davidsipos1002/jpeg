#pragma once

#include <stdio.h>
#include <stdlib.h>

#ifdef __aarch64__
#include <arm_neon.h>
#endif

#define PACKED __attribute__((packed))
#define INLINE __attribute__((always_inline))
#define ENDIAN_SWAP(x) (__builtin_bswap16(x))
#define TO_UINT_ADVANCE(t, x, y) (x) = *((t *) (y)); (y) += sizeof(t)

#define safeMalloc(p, x)           \
do                                 \
{                                  \
    (p) = malloc((x));             \
    if (!p) {                      \
        printf("out of memory\n"); \
        exit(1);                   \
    }                              \
} while(0)

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif

BEGIN_C_DECLS

#ifdef __aarch64__
    typedef float16_t jpegf;
#else
    typedef float jpegf;
#endif

static inline void print_binary16(uint16_t x)
{
    uint16_t mask = 0x8000;
    while(mask)
    {
        printf("%u", (x & mask) != 0);
        mask >>= 1;
    }
}

END_C_DECLS
