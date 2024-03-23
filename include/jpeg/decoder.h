#pragma once

#include <misc.h>

BEGIN_C_DECLS

typedef void* decoder;

#define DUMP_COMPONENT_PARAMS

decoder *init_decoder(const char *filename);
uint8_t decode_image(decoder *dec);
void free_decoder(decoder *dec);

END_C_DECLS
