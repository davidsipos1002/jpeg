#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

void gen_inverse_zigzag();
void dequantize(int16_t *b, uint8_t *q);
void unzigzag(int16_t *b, float **mat);
void idct(float **mat);
void undo_level_shift(float **mat);

END_C_DECLS
