#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

typedef struct
{
    uint32_t count;
    uint32_t row;
    uint32_t col;
    float ***mat;
} arrmat;

arrmat *alloc_matrices(uint32_t count, uint32_t row, uint32_t col);
void free_matrices(arrmat *p);

END_C_DECLS
