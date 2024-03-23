#pragma once

#include <jpeg/format.h>
#include <misc.h>
#include <ds/uthash.h>

BEGIN_C_DECLS

typedef struct
{
    uint8_t class;
    uint8_t destination;
    // encoder tables
    uint8_t ehufsi[256];
    uint16_t ehufco[256];
    // decoder tables
    uint16_t valcount;
    uint8_t *huffval;
    uint16_t mincode[16];
    uint16_t maxcode[16];
    uint16_t valptr[16];
} hfft;

hfft *extract_huffman_table(dht table);
void free_huffman_table(hfft *t);

END_C_DECLS
