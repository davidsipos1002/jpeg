#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

void gen_inverse_zigzag();
void dequantize(int16_t *b, uint8_t *q);
void unzigzag(int16_t *b, float **mat);
void idct(float **mat);
void undo_level_shift(float **mat);
void convert_to_rgb(float y, float cb, float cr, uint8_t *rp, uint8_t *gp, uint8_t *bp);

END_C_DECLS
