#pragma once

#include <misc.h>
#include <jpeg/image.h>

BEGIN_C_DECLS

typedef void* decoder;

decoder *init_decoder(const char *filename);
uint8_t decode_image(decoder *dec);
image *get_image(decoder *dec);
void free_decoder(decoder *dec);

END_C_DECLS