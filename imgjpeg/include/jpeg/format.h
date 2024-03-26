#pragma once

#include <stdint.h>

#include <misc.h>

#define DENSITY_NONE 0x00
#define DENSITY_DPI 0x01
#define DENSITY_DPCM 0x02

#define MRK(x) (0xFF00 | (x))
#define BDCT_SOF 0xC0
#define ESDCT_SOF 0xC1
#define EPDCT_SOF 0xC2
#define DHT 0xC4
#define RST(m) (0xD0 | (m))
#define SOI 0xD8
#define EOI 0xD9
#define SOS 0xDA
#define DQT 0xDB
#define DNL 0xDC
#define DRI 0xDD
#define APP(n) (0xE0 | (n))

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

typedef struct _fhdr_comp
{
    uint8_t c; // component id
    uint8_t h; // horizontal sampling factor
    uint8_t v; // vertical sampling factor
    uint8_t tq; // quantization table selector
    struct _fhdr_comp *next;
    struct _fhdr_comp *prev;
} fhdr_comp;

typedef struct
{
    uint8_t p; // precision
    uint16_t y; // max line count
    uint16_t x; // max sample count per line
    uint8_t nf; // number of components
    fhdr_comp *comp; // component-specification parameters
} fhdr;

typedef struct _shdr_comp
{
    uint8_t cs; // scan component selector
    uint8_t td; // dc entropy coding table destination selector
    uint8_t ta; // ac entropy coding table destination selector
    struct _shdr_comp *next;
    struct _shdr_comp *prev;
} shdr_comp;

typedef struct
{
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
    int16_t *q; // quantization table
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

END_C_DECLS
