#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

void gen_inverse_zigzag();
void dequantize(int16_t *b, int16_t *q);
void unzigzag(int16_t *b, jpegf **mat);
void idct(jpegf **mat);
void undo_level_shift(jpegf **mat);
void convert_to_rgb(jpegf y, jpegf cb, jpegf cr, uint8_t *rp, uint8_t *gp, uint8_t *bp);

END_C_DECLS
