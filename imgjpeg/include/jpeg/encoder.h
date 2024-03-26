#pragma once

#include <misc.h>

#include <jpeg/image.h>

#define JPEG_QUALITY_100 0
#define JPEG_QUALITY_80 1
#define JPEG_QUALITY_50 2

BEGIN_C_DECLS

typedef void* encoder;

encoder *init_encoder(const char *filename, image *img, uint8_t quality);
uint8_t encode_image(encoder *enc);
void free_encoder(encoder *enc);

END_C_DECLS
