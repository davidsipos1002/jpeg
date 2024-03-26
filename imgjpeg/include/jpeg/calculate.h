#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

void gen_inverse_zigzag();
void gen_cos_table();
void dequantize(int16_t *b, int16_t *q);
void unzigzag(int16_t *b, jpegf **mat);
void idct(jpegf **mat);
void undo_level_shift(jpegf **mat);
void convert_to_rgb(jpegf y, jpegf cb, jpegf cr, uint8_t *rp, uint8_t *gp, uint8_t *bp);

void convert_to_ycbcr(jpegf r, jpegf g, jpegf b, jpegf *yp, jpegf *cbp, jpegf *crp);
void fdct(jpegf **mat);
void zigzag(jpegf **mat, int16_t *b);
void quantize(int16_t *b, int16_t *q);

END_C_DECLS
