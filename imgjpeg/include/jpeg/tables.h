#pragma once

#include <misc.h>
#include <jpeg/huffman.h>

BEGIN_C_DECLS

typedef struct
{
    hfft *huffman_dc[4];
    hfft *huffman_ac[4];
    int16_t *quantization[4];
} tables;

tables *init_tables();
void install_huffman_table(tables *t, hfft *h);
void install_quantization_table(tables *t, dqt *q);
void free_tables(tables *t);

END_C_DECLS
