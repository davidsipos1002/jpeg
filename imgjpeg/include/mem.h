#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

typedef struct
{
    uint32_t count;
    uint32_t row;
    uint32_t col;
    jpegf ***mat;
} arrmat;

arrmat *alloc_matrices(uint32_t count, uint32_t row, uint32_t col);
void free_matrices(arrmat *p);

uint8_t **alloc_mat(uint16_t row, uint16_t col);
void free_mat(uint8_t **mat, uint16_t row, uint16_t col);

END_C_DECLS
