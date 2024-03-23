#pragma once

#include <stdint.h>

#include <misc.h>

#define DENSITY_NONE 0x00
#define DENSITY_DPI 0x01
#define DENSITY_DPCM 0x02

BEGIN_C_DECLS

typedef struct
{
    uint16_t lp; // total APP0 count
    char identifier[5]; // "JFIF\0"
    uint16_t version;
    uint8_t units; // denstity units 
    uint16_t hdensity; // horizontal density
    uint16_t vdensity; // vertical density
    uint16_t hthumbnail; // horizontal thumbnail pixel count
    uint16_t vthumbnail; // vertical thumbnail pixel count
    uint8_t rgb[];
} jfif_app0;

typedef struct
{
    uint8_t c; // component id
    uint8_t h; // horizontal sampling factor
    uint8_t v; // vertical sampling factor
    uint8_t tq; // quantization table selector
} fhdr_comp;

typedef struct
{
    uint16_t lf; // frame header length
    uint8_t p; // precision
    uint16_t y; // max line count
    uint16_t x; // max sample count per line
    uint8_t nf; // number of components
    fhdr_comp *comp; // component-specification parameters
} fhdr;

typedef struct
{
    uint8_t cs; // scan component selector
    uint8_t td; // dc entropy coding table destination selector
    uint8_t ta; // ac entropy coding table destination selector
} shdr_comp;

typedef struct
{
    uint16_t ls; // scan header length
    uint8_t ns; // number of components in the scan
    shdr_comp *comp; // component-specification parameters
    uint8_t ss; // start of spectral selection
    uint8_t se; // end of spectral selection
    uint8_t ah; // successive approximation bit position high
    uint8_t al; // successive approximation bit position low
} shdr;

#define QUANTIZATION_8 0
#define QUANTIZATION_16 1
typedef struct _dqt
{
    uint8_t pq; // quantization table element precision
    uint8_t tq; // quantization table destination identifier
    uint8_t *q; // quantization table
    struct _dqt *next;
    struct _dqt *prev;
} dqt;

#define HUFFMAN_DC 0
#define HUFFMAN_AC 1

typedef struct _dhtt
{
    uint8_t tc; // table class
    uint8_t th; // Huffman table destination identifier
    uint8_t l[16]; // number of Huffman codes of length i
    uint16_t lv; // length of v
    uint8_t *v; // values associated with each Huffman code
    struct _dhtt *next;
    struct _dhtt *prev;
} dhtt;

typedef dhtt* dht;
typedef struct 
{
    uint16_t lr; // define length of restart segment length, always 4
    uint16_t ri; // number of MCUs in a restart interval
} dri;

typedef struct 
{
    uint16_t ld; // define number of lines segment length, always 4
    uint16_t nl; // number of lines
} dnl;

END_C_DECLS