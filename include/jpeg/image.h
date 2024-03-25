#pragma once

#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

typedef struct
{
    uint16_t y;
    uint16_t x;
    uint8_t **r;
    uint8_t **g;
    uint8_t **b;
} image;

image *init_image(uint16_t rows, uint16_t cols);
void free_image(image *img);

END_C_DECLS
