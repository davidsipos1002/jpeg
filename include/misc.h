#pragma once

#define PACKED __attribute__((packed))
#define INLINE __attribute__((always_inline))
#define ENDIAN_SWAP(x) (__builtin_bswap16(x))