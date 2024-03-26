#include <jpeg/encoder.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <jpeg/calculate.h>
#include <jpeg/format.h>
#include <jpeg/huffman.h>

#include <mem.h>

typedef struct
{
    FILE *f; // file to which we are writing
    image *img; // image to encode
    uint8_t q; // quality
    uint16_t py; // padded number of lines (multiple of 8)
    uint16_t px; // padded number of columns (multiple of 8)
    uint16_t blcy; // number of block-rows
    uint16_t blcx; // number of blocks on a row 
    uint32_t blc; // total number of blocks in a component
    arrmat *bl[3]; // image converted to data units
    int16_t **datau[3]; // data units to entropy encode  
    hfft *tbl[3][2]; // huffman tables used for encoding
    int16_t pred[3]; // predictions for DC
    uint8_t b; // byte waiting to be written
    uint8_t cnt; // used to where are we in the current byte
} encoder_s;

// very simple encoder, meaning fixed Huffman tables, no chroma subsampling, no restart
// intervals and only 3 quality factors 100, 80, 50

// quantization tables taken from images in the img directory
static int16_t qtables[3][2][64] = 
{
    // quality 100
    {
        // Luma
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
        // Chroma
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
    },
    // quality 80
    {
        // Luma
        {6, 4, 5, 6, 5, 4, 6, 6, 5, 6, 7, 7, 6, 8, 10, 16, 10, 
        10, 9, 9, 10, 20, 14, 15, 12, 16, 23, 20, 24, 24, 23, 
        20, 22, 22, 26, 29, 37, 31, 26, 27, 35, 28, 22, 22, 32, 
        44, 32, 35, 38, 39, 41, 42, 41, 25, 31, 45, 
        48, 45, 40, 48, 37, 40, 41, 40},
        // Chroma
        {7, 7, 7, 10, 8, 10, 19, 10, 10, 19, 40, 26, 22, 26, 40,
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
         40, 40, 40, 40, 40, 40, 40}
    },
    // quality 50, they can also be found in the standard
    {
        // Luma
        {16, 11, 12, 14, 12, 10, 16, 14, 13, 14, 18, 17, 16, 19, 24,
        40, 26, 24, 22, 22, 24, 49, 35, 37, 29, 40, 58, 51, 61, 60,
        57, 51, 56, 55, 64, 72, 92, 78, 64, 68, 87, 69, 55, 56, 80,
        109, 81, 87, 95, 98, 103, 104, 103, 62, 77, 113, 121, 112,
        100, 120, 92, 101, 103, 99},
        // Chroma
        {17, 18, 18, 24, 21, 24, 47, 26, 26, 47, 99, 66, 56, 66, 99,
        99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
        99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 
        99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99,
        99, 99, 99, 99, 99}
    }
};

// Huffman tables taken from the specification

// Luma tables
static uint8_t lumadcbits[] =
{
    0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0
};

static uint8_t lumadchuffval[] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
};

static uint8_t lumaacbits[] =
{
    0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 125
};

static uint8_t lumaachuffval[] = 
{
    1, 2, 3, 0, 4, 17, 5, 18, 33, 49, 65, 6, 19, 81, 97,
    7, 34, 113, 20, 50, 129, 145, 161, 8, 35, 66, 177, 193,
    21, 82, 209, 240, 36, 51, 98, 114, 130, 9, 10, 22, 23,
    24, 25, 26, 37, 38, 39, 40, 41, 42, 52, 53, 54, 55, 56,
    57, 58, 67, 68, 69, 70, 71, 72, 73, 74, 83, 84, 85, 86,
    87, 88, 89, 90, 99, 100, 101, 102, 103, 104, 105, 106, 
    115, 116, 117, 118, 119, 120, 121, 122, 131, 132, 133, 
    134, 135, 136, 137, 138, 146, 147, 148, 149, 150, 151, 
    152, 153, 154, 162, 163, 164, 165, 166, 167, 168, 169, 
    170, 178, 179, 180, 181, 182, 183, 184, 185, 186, 194, 
    195, 196, 197, 198, 199, 200, 201, 202, 210, 211, 212, 
    213, 214, 215, 216, 217, 218, 225, 226, 227, 228, 229, 
    230, 231, 232, 233, 234, 241, 242, 243, 244, 245, 246, 
    247, 248, 249, 250
};

// Chroma tables
static uint8_t chromadcbits[] = 
{
    0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 
};

static uint8_t chromadchuffval[] =
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
};

static uint8_t chromaacbits[] =
{
    0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 119
};

static uint8_t chromaachuffval[] =
{
    0, 1, 2, 3, 17, 4, 5, 33, 49,
    6, 18, 65, 81, 7, 97, 113, 19, 
    34, 50, 129, 8, 20, 66, 145, 161, 177, 
    193, 9, 35, 51, 82, 240, 21, 98, 114, 
    209, 10, 22, 36, 52, 225, 37, 241, 23, 
    24, 25, 26, 38, 39, 40, 41, 42, 53, 
    54, 55, 56, 57, 58, 67, 68, 69, 70, 71, 
    72, 73, 74, 83, 84, 85, 86, 87, 88, 89,
    90, 99, 100, 101, 102, 103, 104, 105,
    106, 115, 116, 117, 118, 119, 120, 121, 
    122, 130, 131, 132, 133, 134, 135, 136, 
    137, 138, 146, 147, 148, 149, 150, 151, 
    152, 153, 154, 162, 163, 164, 165, 166, 167, 
    168, 169, 170, 178, 179, 180, 181, 182, 183, 
    184, 185, 186, 194, 195, 196, 197, 198, 199, 
    200, 201, 202, 210, 211, 212, 213, 214, 215, 216, 
    217, 218, 226, 227, 228, 229, 230, 231, 232, 233, 
    234, 242, 243, 244, 245, 246, 247, 248, 249, 250
};

encoder *init_encoder(const char *filename, image *img, uint8_t quality)
{
    encoder_s *ret;
    safeMalloc(ret, sizeof(encoder_s));
    memset(ret, 0, sizeof(encoder_s));
    ret->f = fopen(filename, "wb");
    ret->img = img;
    ret->q = quality;
    return (encoder *) ret;
}

static void write_dqt(encoder_s *e)
{
    static uint8_t dqt_header[] = 
    {
        0xFF, 0xDB, // DQT marker
        0x00, 0x43, // segment length
        0x00 // precision and destination
    };
    
    // write Luma qt
    dqt_header[4] = 0x00;
    fwrite(dqt_header, sizeof(dqt_header), 1, e->f);
    int16_t *tbl = qtables[e->q][0];
    uint8_t x;
    for (uint8_t i = 0; i < 64;i++)
    {
        x = (uint8_t) tbl[i];
        fwrite(&x, 1, 1, e->f);
    }
    
    // write Chroma qt
    dqt_header[4] = 0x01;
    fwrite(dqt_header, sizeof(dqt_header), 1, e->f);
    tbl = qtables[e->q][1];
    for (uint8_t i = 0; i < 64;i++)
    {
        x = (uint8_t) tbl[i];
        fwrite(&x, 1, 1, e->f);
    }
    
}

static void write_sof(encoder_s *e)
{
    static uint8_t frame_hdr[] =
    {
        0xFF, 0xC0, // SOF marker for BDCT
        0x00, 0x11, // lf
        0x08, // precision 
        0x00, 0x00, // Y, number lines
        0x00, 0x00, // X, number of samples per line
        0x03, // Nf, number of image components
        0x01, 0x11, 0x00, // Luma, component 1, no subsampling and use qt 0
        0x02, 0x11, 0x01, // Cb, component 2, no subsampling and use qt 1
        0x03, 0x11, 0x01 // Cr, component 3, nu subsampling and use qt1
    };

    // fill image size
    frame_hdr[5] = e->img->y >> 8;
    frame_hdr[6] = e->img->y & 0xFF;
    frame_hdr[7] = e->img->x >> 8;
    frame_hdr[8] = e->img->x & 0xFF;
    // write frame header
    fwrite(frame_hdr, sizeof(frame_hdr), 1, e->f);
}

static void write_dht(encoder_s *e)
{
    static uint8_t huffman_hdr[] =
    {
        0xFF, 0xC4, // DHT marker
        0x00, 0x00, // lh length
        0x00, // tc, table class 0 = DC, 1 = AC and tq, destination 0-3
    };
    
    // write Luma huffman tables

    // write DC tc=0, th=0
    uint16_t lh = 19 + sizeof(lumadchuffval);
    huffman_hdr[2] = lh >> 8;
    huffman_hdr[3] = lh & 0xFF;
    huffman_hdr[4] = 0;
    fwrite(huffman_hdr, sizeof(huffman_hdr), 1, e->f);
    fwrite(lumadcbits, sizeof(lumadcbits), 1, e->f);
    fwrite(lumadchuffval, sizeof(lumadchuffval), 1, e->f);
    
    // write AC tc=1 th=0
    lh = 19 + sizeof(lumaachuffval);
    huffman_hdr[2] = lh >> 8;
    huffman_hdr[3] = lh & 0xFF;
    huffman_hdr[4] = 0x10;    
    fwrite(huffman_hdr, sizeof(huffman_hdr), 1, e->f);
    fwrite(lumaacbits, sizeof(lumaacbits), 1, e->f);
    fwrite(lumaachuffval, sizeof(lumaachuffval), 1, e->f);
    
    // write Chroma huffman tables
    
    // write DC tc=0 th=1
    lh = 19 + sizeof(chromadchuffval);
    huffman_hdr[2] = lh >> 8;
    huffman_hdr[3] = lh & 0xFF;
    huffman_hdr[4] = 0x01;
    fwrite(huffman_hdr, sizeof(huffman_hdr), 1, e->f);
    fwrite(chromadcbits, sizeof(chromadcbits), 1, e->f);
    fwrite(chromadchuffval, sizeof(chromadchuffval), 1, e->f);
    
    // write AC tc=1 th=1
    lh = 19 + sizeof(chromaachuffval);
    huffman_hdr[2] = lh >> 8;
    huffman_hdr[3] = lh & 0xFF;
    huffman_hdr[4] = 0x11;    
    fwrite(huffman_hdr, sizeof(huffman_hdr), 1, e->f);
    fwrite(chromaacbits, sizeof(chromaacbits), 1, e->f);
    fwrite(chromaachuffval, sizeof(chromaachuffval), 1, e->f);
}

static void write_sos(encoder_s *e)
{
    static uint8_t scan_header[] =
    {
        0xFF, 0xDA, // SOS marker
        0x00, 0x0C, // header length
        0x03, // Ns, number of components in the scan, we do interleaved scan
        0x01, 0x00, // Luma, component 1, td=0 ta=0
        0x02, 0x11, // Cb, component 2, td=1 ta=1
        0x03, 0x11, // Cr, component 3, td=1 ta=1
        0x00, // spectral start we don't do progressive but we set them anyways
        0x3F, // spectral end
        0x00, // ah=0, al=0 
    };
    
    // write the scan header
    fwrite(scan_header, sizeof(scan_header), 1, e->f);
}

static void pad_component(encoder_s *e, uint8_t **comp, uint8_t c)
{
    arrmat *m = alloc_matrices(e->blc, 8, 8);
    uint32_t bndx = 0;
    uint32_t i, j;
    for (uint16_t x = 0; x < e->blcy; x++)
    {
        for (uint16_t y = 0; y < e->blcx; y++)
        {
            i = x * 8;
            for (uint8_t k = 0; k < 8; k++)
            {
                j = y * 8;
                for (uint8_t l = 0; l < 8; l++)
                {
                    m->mat[bndx][k][l] = (jpegf) comp[i][j];
                    if (j != e->img->x - 1)
                        j++;
                }
                if (i != e->img->y - 1)
                    i++;
            }
            bndx++;
        }
    }
    e->bl[c] = m;
}

static void convert_components(encoder_s *e)
{
    for (uint32_t b = 0; b < e->blc; b++)
    {
        // TODO: this can be optimized
        jpegf **rb = e->bl[0]->mat[b];
        jpegf **gb = e->bl[1]->mat[b];
        jpegf **bb = e->bl[2]->mat[b];   
        for (uint8_t i = 0; i < 8; i++)
        {
            for (uint8_t j = 0; j < 8; j++)
            {
                jpegf y, cb, cr;
                convert_to_ycbcr(rb[i][j], gb[i][j], bb[i][j], &y, &cb, &cr);
                // do the level shift
                y -= 128;
                cb -= 128;
                cr -= 128;
                rb[i][j] = y;
                gb[i][j] = cb;
                bb[i][j] = cr;
            }
        }
    }
}

static void prepare_component(encoder_s *e, uint8_t comp, uint8_t chroma)
{
    safeMalloc(e->datau[comp], e->blc * sizeof(int16_t *));
    for (uint32_t b = 0; b < e->blc; b++)
    {
        safeMalloc(e->datau[comp][b], 64 * sizeof(int16_t));
        fdct(e->bl[comp]->mat[b]);
        zigzag(e->bl[comp]->mat[b], e->datau[comp][b]);
        quantize(e->datau[comp][b], qtables[e->q][chroma]);
    }
    free_matrices(e->bl[comp]);
}

static void prepare_for_entropy_coding(encoder_s *e)
{
    prepare_component(e, 0, 0);
    prepare_component(e, 1, 1);
    prepare_component(e, 2, 1);
}

static void load_huffman(encoder_s *e)
{
    dhtt currtbl;
    currtbl.next = NULL;
    currtbl.prev = NULL;
    hfft *hft;
    // load Luma tables 

    printf("loading huffman in the encoder\n"); 
    // Luma DC tc = 0 th = 0
    currtbl.tc = 0;
    currtbl.th = 0;
    memcpy(currtbl.l, lumadcbits, 16 * sizeof(uint8_t));
    currtbl.lv = sizeof(lumadchuffval);
    currtbl.v = lumadchuffval;
    hft = extract_huffman_table(&currtbl);
    e->tbl[0][0] = hft;

    // Luma AC tc = 1 th = 0 
    currtbl.tc = 1;
    currtbl.th = 0;
    memcpy(currtbl.l, lumaacbits, 16 * sizeof(uint8_t));    
    currtbl.lv = sizeof(lumaachuffval);
    currtbl.v = lumaachuffval;
    hft = extract_huffman_table(&currtbl);
    e->tbl[0][1] = hft;

    // Chroma DC tc = 0 th = 1
    currtbl.tc = 0;
    currtbl.th = 1;
    memcpy(currtbl.l, chromadcbits, 16 * sizeof(uint8_t));    
    currtbl.lv = sizeof(chromadchuffval);
    currtbl.v = chromadchuffval;
    hft = extract_huffman_table(&currtbl);
    e->tbl[1][0] = hft;
    e->tbl[2][0] = hft;
    
    // Chroma AC tc = 1 th = 1
    currtbl.tc = 1;
    currtbl.th = 1;
    memcpy(currtbl.l, chromaacbits, 16 * sizeof(uint8_t));    
    currtbl.lv = sizeof(chromaachuffval);
    currtbl.v = chromaachuffval;
    hft = extract_huffman_table(&currtbl);
    e->tbl[1][1] = hft;
    e->tbl[2][1] = hft;
}

static void writebit(encoder_s *e, uint8_t bit)
{
    static const uint8_t zero = 0x00;
    if (e->cnt == 8)
    {
        fwrite(&(e->b), 1, 1, e->f);
        if (e->b == 0xFF)
            fwrite(&zero, 1, 1, e->f);
        e->cnt = 0;
        e->b = 0;
    }
    e->b = (e->b << 1) | (bit != 0);
    e->cnt++;
}

static uint8_t get_cat(int16_t x, int16_t mxcat)
{
    if (!x)
        return 0;
    int16_t absx = x < 0 ? -x : x;
    for (int16_t cat = 1; cat <= mxcat; cat++)
    {
        int16_t low = ((int16_t) 1 << (cat - 1)); 
        int16_t high = ((int16_t) 1 << (cat)) - 1;
        if (low <= absx && absx <= high)
            return cat;
    }
    return mxcat;
}

static void write_bits(encoder_s *e, uint16_t x, uint8_t count)
{
    // this not expected by the function
    assert(count != 0);
    uint16_t mask = (uint16_t) 1 << (count - 1); 
    while (mask)
    {
        uint8_t bit = (x & mask) != 0;
        writebit(e, bit);
        mask >>= 1;
    }
}

static void encode_dc(encoder_s *e, uint8_t comp, int16_t diff)
{
    hfft *h = e->tbl[comp][0];
    uint8_t cat = get_cat(diff, 11);
    uint8_t sz = h->ehufsi[cat];
    uint16_t code = h->ehufco[cat];
    write_bits(e, code, sz);
    if (diff)
    {
        if (diff < 0)
            diff--;
        uint16_t x = *((uint16_t *) &diff); 
        write_bits(e, x, cat);
    }
}

static void encode_ac(encoder_s *e, uint8_t comp, uint8_t r, int16_t ac)
{
    hfft *h = e->tbl[comp][1];
    uint8_t cat = get_cat(ac, 10); 
    uint8_t rs = (r << 4) | (cat & 0xF);
    uint8_t sz = h->ehufsi[rs];
    uint16_t code = h->ehufco[rs]; 
    write_bits(e, code, sz);
    if (ac < 0)
        ac--;
    uint16_t x = *((uint16_t *) &ac);
    write_bits(e, x, cat);
}

static void write_special(encoder_s *e, uint8_t comp, uint8_t val)
{
    assert(val == 0x00 || val == 0xF0);
    hfft *h = e->tbl[comp][1];
    uint8_t sz = h->ehufsi[val];
    uint16_t code = h->ehufco[val];
    write_bits(e, code, sz);
}

static void entropy_encode(encoder_s *e)
{
    for (uint32_t b = 0; b < e->blc; b++)
    {
        for (uint8_t comp = 0; comp < 3; comp++)
        {   
            // encode DC coefficient
            int16_t *block = e->datau[comp][b];
            int16_t diff = block[0] - e->pred[comp];
            encode_dc(e, comp, diff);
            e->pred[comp] = block[0];
            
            // encode AC coefficients
            uint8_t r = 0;
            for (uint8_t k = 1; k < 64; k++)
            {
                if (!block[k])
                {
                    // write EOB
                    if (k == 63)
                        write_special(e, comp, 0x00);
                    else
                        r++;
                    continue;
                }
                // write ZRL
                while (r > 15)
                {
                    write_special(e, comp, 0xF0);
                    r -= 16;
                }
                encode_ac(e, comp, r, block[k]);
                r = 0;
            }
        } 
    }
    // we might remain with a complete byte which is not written
    // we need another call to flush it
    // or if we have an incomplete byte we need to complete it
    // so that the EOI will be at byte boundary
    if (e->cnt && e->cnt < 8)
    {
        uint8_t pad = 8 - e->cnt;
        for (uint8_t i = 0; i < pad; i++)
            writebit(e, 1); 
    }
    if (e->cnt == 8)
        writebit(e, 0);
}

uint8_t encode_image(encoder *enc)
{
    static uint8_t encoder_app0[] = 
    { 
        0xFF, 0xE0, // APP0 marker
        0x00, 0x10, // lp
        0x4A, 0x46, 0x49, 0x46, 0x00, // identifier "JFIF\0"
        0x01, 0x01, // version
        0x00, // unspecified units meaning aspect ratio y/x = vdensity/hdensity
        0x00, 0x00, // h density
        0x00, 0x00, // v density
        0x00, // HthumbnailA, no thumbnail
        0x00 // VthumbnailA, no thumbnail
    };

    encoder_s *e = (encoder_s *) enc;

    // write start of image marker
    uint16_t marker = MRK(SOI);
    marker = ENDIAN_SWAP(marker);
    fwrite(&marker, sizeof(uint16_t), 1, e->f); 

    // write APP0 marker
    encoder_app0[12] = e->img->x >> 8;
    encoder_app0[13] = e->img->x & 0xFF;
    encoder_app0[14] = e->img->y >> 8;
    encoder_app0[15] = e->img->y & 0xFF;
    fwrite(&encoder_app0, sizeof(encoder_app0), 1, e->f); 
    
    // write quantization tables
    write_dqt(e);

    // write frame header
    write_sof(e);

    // write huffman tables
    write_dht(e); 
    
    // write scan header
    write_sos(e);

    // steps to perform:
    // 1. pad the image so that we have dimensions multiples of 8
    // and convert components to blocks
    // 2. convert to YCbCr
    // 3. FDCT
    // 4. Zigzag
    // 5. Quantize
    // 6. Entropy encode  
    
    uint16_t py = e->img->y;
    if (py % 8)
        py += 8 - (py % 8);
    e->py = py;
    e->blcy = py / 8;
    uint16_t px = e->img->x;
    if (px % 8)
        px += 8 - (px % 8);
    e->px = px;
    e->blcx = px / 8;
    e->blc = e->blcy * e->blcx;
    
    pad_component(e, e->img->r, 0);
    pad_component(e, e->img->g, 1);
    pad_component(e, e->img->b, 2);
    
    convert_components(e); 
    
    prepare_for_entropy_coding(e);
    
    load_huffman(e);

    entropy_encode(e);
    
    // write end of image marker
    marker = MRK(EOI);
    marker = ENDIAN_SWAP(marker);
    fwrite(&marker, sizeof(uint16_t), 1, e->f); 
    
    return 1;
}

void free_encoder(encoder *enc)
{
    encoder_s *e = (encoder_s *) enc;
    fclose(e->f);
    free_huffman_table(e->tbl[0][0]);
    free_huffman_table(e->tbl[0][1]);
    free_huffman_table(e->tbl[1][0]);
    free_huffman_table(e->tbl[1][1]);
    for (uint8_t i = 0; i < 3; i++)
    {
        for (uint32_t b = 0; b < e->blc; b++)
            free(e->datau[i][b]);
        free(e->datau[i]);
    }
    free(e);
}
