#include <jpeg/decoder.h>

#include <assert.h>
#include <math.h>

#include <file.h>
#include <mem.h>
#include <threads.h>
#include <jpeg/tables.h>
#include <jpeg/parser.h>
#include <jpeg/calculate.h>
#include <jpeg/image.h>

#define MODE_BDCT 0

typedef struct
{
    uint8_t h; // horizontal sampling factor
    uint8_t v; // vertical sampling factor
    uint8_t tq; // quantization table selector
    uint16_t y; // number of lines
    uint16_t x; // number of columns
    uint16_t py; // padded number of lines (multiple of 8)
    uint16_t px; // padded number of columns (multiple of 8)
    uint16_t blcy; // number of block-rows (multiple of v)
    uint16_t blcx; // number of blocks on a row (multiple of h)
    uint32_t mcuc; // number of mcus according to sampling factors
    uint32_t blc; // total block count
    int16_t **blocks; // decoded blocks of the component
    arrmat *m; // array of 8x8 matrices
} component;

typedef struct
{
    uint8_t ns; // number of components in the scan
    uint8_t comps[4]; // component in the current scan
    uint8_t tbl[3][2]; // huffman tables used for decoding compoents in the scan
    uint8_t ss; // spectral start
    uint8_t se; // spectral end
    uint8_t ah; // successive approx bit pos high
    uint8_t al; // successive approx bit pos low
    uint16_t ricount; // number of restart intervals in the scan
    uint16_t ril; // mcus in a restart interval
    uint32_t mcur; // remaining mcus to decode
    uint8_t currc; // current component
    uint32_t mcucurr; // current mcu 
    int16_t pred[3]; // DC prediction
    uint8_t cnt; // used by nextbit
} scan;

typedef struct
{
    mmapfile *file; // handle to memory-file
    tables *tables; // installed decoder tables
    uint16_t ril; // restart interval length
    uint8_t *p; // current file position
    uint8_t mode; // decoder mode
    uint16_t y; // max number of lines
    uint16_t x; // max number of columns
    uint8_t nf; // number of image components
    uint16_t hmax; // max horizontal sampling factor
    uint16_t vmax; // max vertical sampling factor
    component comps[3]; // component information
    scan currscan; // current scan information
    uint8_t gotimg; // remembers whether the image pointer was requested
    image *img; // final image
} decoder_s;

// thread parameter 
typedef struct
{
    decoder_s *d; // decoder instance
    uint8_t comp; // component to process
    pthread_mutex_t *m; // mutex protecting rdy
    pthread_cond_t *c; // condition variable associated with the mutex
    volatile uint8_t *rdy; // how many threads are finished
} decoder_tp;

// calculate component properties which directly the decoding process
static void compute_component_params(decoder_s *d, component *c)
{
    // calculate the component's dimension
    uint16_t comp_y = ceil(d->y * (double) c->v / d->vmax);
    c->y = comp_y;
    uint16_t comp_x = ceil(d->x * (double) c->h / d->hmax);
    c->x = comp_x;
    // padding to make an integer number of data units
    c->py = c->y;
    if (c->y % 8)
        c->py += 8 - (c->y % 8);
    c->px = c->x;
    if (c->x % 8)
        c->px += 8 - (c->x % 8);
    // pad a row to make the data unit count multiple of h 
    uint16_t blx = c->px / 8;
    c->blcx = blx;
    if (blx % c->h)
        c->blcx += c->h - (blx % c->h);
    // pad the component to make an integer number of block-rows
    uint16_t bly = c->py / 8;
    c->blcy = bly;
    if (bly % c->v)
        c->blcy += c->v - (bly % c->v);
    // allocate memory for the blocks of the component
    uint32_t block_count = c->blcx * c->blcy; 
    safeMalloc(c->blocks, block_count * sizeof(int16_t *));
    for (uint32_t i = 0; i < block_count; i++)
    {
        safeMalloc(c->blocks[i], 64 * sizeof(int16_t));
        memset(c->blocks[i], 0, 64 * sizeof(int16_t));
    }
    c->blc = block_count;
    c->mcuc = (c->blcy / c->v) * (c->blcx / c->h);
}

// allocate memory for the decoder
decoder *init_decoder(const char *filename)
{
    decoder_s *ret;
    safeMalloc(ret, sizeof(decoder_s));
    memset(ret, 0, sizeof(decoder_s));
    ret->file = map_file(filename);
    ret->p = ret->file->p;
    ret->tables = init_tables();
    // return the decoder struct as opaque pointer to hide the struct layout
    return (decoder *) ret;
}

// reads two bytes and interprets them as a marker
static uint16_t get_marker(decoder_s *d)
{
    uint16_t marker;
    TO_UINT_ADVANCE(uint16_t, marker, d->p); 
    return ENDIAN_SWAP(marker);
}

// interpret markers: DHT (define huffman table), DQT (define quantization table)
// DRI (define restart interval)
// for any other marker jump over it
static uint8_t interpret_markers(decoder_s *d, uint16_t marker)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Interpreting marker %x\n", marker);
#endif
    uint16_t len;
    uint8_t stat = 1;
    switch (marker)
    {
        case MRK(DHT): 
        {
            dht huff = parse_dht(d->p, &len); 
            dht el;
            DL_FOREACH(huff, el)
            {
                hfft *table = extract_huffman_table(el);
                install_huffman_table(d->tables, table); 
            }
            free_dht(huff);
            break;
        }
        case MRK(DQT):
        {
            dqt *quant = parse_dqt(d->p, &len);
            dqt *curr;
            DL_FOREACH(quant, curr)
                install_quantization_table(d->tables, curr);
            free_dqt(quant);
            break;
        } 
        case MRK(DRI):
        {
            d->ril = parse_dri(d->p, &len);
            break;
        }
        case MRK(APP(0)):
        {
            stat = parse_app0(d->p, &len);
            break;
        }
        default:
        {
            len = get_marker_seg_len(d->p);
            break;
        }
    }
    d->p += len;
    return stat;
}

// extend the sign bit of a decoded value V 
static int16_t extend(int16_t v, uint8_t t)
{
    if (!t)
        return 0;
    int16_t vt = 1 << (t - 1); 
    if (v < vt)
    {
        vt = -1;
        vt <<= t;
        vt++;
        return v + vt;
    }
    return v;
}

// get the next bit of the file, ignore byte stuffing
static uint8_t nextbit(decoder_s *d)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Nextbit\n");
#endif
    static uint8_t b;
    if (!d->currscan.cnt)
    {
        TO_UINT_ADVANCE(uint8_t, b, d->p);
#ifdef DECODER_LOG
        printf("[Decoder]: Read new byte %X from %llX\n", b,(uint64_t) ((d->p - sizeof(uint8_t)) - (uint8_t *) d->file->p));
#endif
        d->currscan.cnt = 8;
        if (b == 0xFF)
        {
            uint8_t b2 = TO_UINT_ADVANCE(uint8_t, b2, d->p); 
            // we don't expect any marker here because we don't support DNL
            assert(b2 == 0);
        }
    }
    uint8_t ret = b >> 7;
    d->currscan.cnt--;
    b <<= 1;
    return ret;
}

// receive s bits from the file
static int16_t receive(decoder_s *d, uint8_t s)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Receiving %u bits from the file\n", s);
#endif
    int16_t v = 0;
    for (uint8_t i = 0;i < s; i++)
        v = (v << 1) + nextbit(d);
    return v;
}

// fetch and interpret the next Huffman code 
uint8_t decode(decoder_s *d, uint8_t coeff)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Interpreting next Huffman code\n");
#endif
    uint8_t i = 0;
    uint16_t code = nextbit(d);
    uint8_t currc = d->currscan.currc;
    uint8_t tsel = d->currscan.tbl[currc][coeff];
    hfft *hft;
    if (coeff == HUFFMAN_DC)
        hft = d->tables->huffman_dc[tsel];
    else
        hft = d->tables->huffman_ac[tsel];
    while (code > hft->maxcode[i] || (hft->maxcode[i] == hft->mincode[i] && hft->maxcode[i] == 0xFFFF))
    {
        i++;
        code = (code << 1) + nextbit(d); 
    }
    uint16_t j = hft->valptr[i];
    j = j + code - hft->mincode[i];
#ifdef DECODER_LOG
    printf("[Decoder]: Decoded => code %u valptr %u size %u val %u", code, hft->valptr[i], i, hft->huffval[j]);
    printf(" huffman code: ");
    print_binary16(code);
    printf("\n");
#endif
    return hft->huffval[j];
}

static uint8_t decode_data_unit(decoder_s *d, int16_t *block)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding new data unit\n");
#endif
    // decode DC DIFF
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding DC...\n");
#endif
    uint8_t t = decode(d, HUFFMAN_DC);
#ifdef DECODER_LOG
    printf("[Decoder]: Decode return size %u\n", t);
#endif
    int16_t diff = receive(d, t);
#ifdef DECODER_LOG
    printf("[Decoder]: Received %d\n", diff);
#endif
    // obtain DC coefficient
    diff = extend(diff, t);
    block[0] = diff + d->currscan.pred[d->currscan.currc];
    d->currscan.pred[d->currscan.currc] = block[0];
    
    // decode AC coefficients
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding AC...\n");
#endif
    uint8_t k = 1;
    while (k < 64)
    {
#ifdef DECODER_LOG
        printf("k is %u\n", k);
#endif
        uint8_t rs = decode(d, HUFFMAN_AC);
        uint8_t s = rs & 0xF;
        uint8_t r = rs >> 4;
#ifdef DECODER_LOG
        printf("[Decoder]: Decode returned %u s %u r %u\n", rs, s, r);
#endif
        if (!s)
        {
            if (r == 15)
            {
                k += 16;
                continue;
            } else
                return 1;
        }
        k += r;
        block[k] = receive(d, s);
#ifdef DECODER_LOG
        printf("[Decoder]: Received %u\n", block[k]);
#endif
        block[k] = extend(block[k], s);
        k++;
    }
    return 1;
}

static uint8_t decode_mcu(decoder_s *d)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding MCU\n");
#endif
    uint8_t cpos = 0;
    uint8_t c;
    while ((c = d->currscan.comps[cpos]) != 0xFF)
    {
        component *comp = &d->comps[c]; 
        uint16_t mcu_rows = comp->blcy / comp->v;
        uint16_t mcu_cols = comp->blcx / comp->h;
        uint16_t mcuy = d->currscan.mcucurr / mcu_cols;
        d->currscan.currc = c;
#ifdef DECODER_LOG
        printf("[Decoder]: Current component %u td %u ta %u\n", c, d->currscan.tbl[c][HUFFMAN_DC], d->currscan.tbl[c][HUFFMAN_AC]);
#endif
        if (mcuy >= mcu_rows)
            return 0;
        uint16_t mcux = d->currscan.mcucurr % mcu_cols; 
        for (uint8_t i = 0; i < comp->v; i++)
        {
            for (uint8_t j = 0; j < comp->h; j++)
            {
                uint16_t bi = mcuy * comp->v + i;
                uint16_t bj = mcux * comp->h + j;
                int16_t *block = comp->blocks[bi * comp->blcx + bj];
                decode_data_unit(d, block);
            }
        }
        cpos++; 
    }
    d->currscan.mcucurr++;
    return 1;
}

// this function checks if its input is a marker, i.e. if it of the form
// 0xFFLL where LL != 0
static uint8_t is_marker(uint16_t marker)
{
    return (marker >> 8) == 0xFF && (marker & 0xFF) != 0;
}

// this decodes a new restart interval
static uint8_t decode_restart_interval(decoder_s *d)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding new restart interval\n");
#endif
    // reset decoder
    memset(d->currscan.pred, 0, 3 * sizeof(int16_t));
    d->currscan.cnt = 0;
    // mcus to decode in this restart interval
    uint32_t mcud = d->currscan.mcur >= d->currscan.ril ? d->currscan.ril : d->currscan.mcur;    
    for (uint32_t i = 0; i < mcud; i++)
    {
        if (!decode_mcu(d))
            return 0;
    } 
    d->currscan.mcur -= mcud;
    uint16_t marker; 
#ifdef DECODER_LOG
    printf("[Decoder]: Restart interval finished we are at %p\n", (uint8_t *) (d->p - (uint8_t *) d->file->p));
#endif
    do
    {
        marker = get_marker(d);
    } while (!is_marker(marker));
    if (marker != MRK(RST(0)) && marker != MRK(RST(1)) && marker != MRK(RST(2)) && 
        marker != MRK(RST(3)) && marker != MRK(RST(4)) && marker != MRK(RST(5)) && 
        marker != MRK(RST(6)) && marker != MRK(RST(7)))
    {
        d->p -= sizeof(uint16_t); 
    }
    return 1;
}

// this function decodes an image scan
static uint8_t decode_scan(decoder_s *d)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding new scan\n");
#endif
    uint16_t len;
    shdr *s = parse_shdr(d->p, &len);
    d->p += len;
    d->currscan.ns = s->ns;
    d->currscan.ss = s->ss;
    d->currscan.se = s->se;
    d->currscan.cnt = 0;
    for (uint8_t i = 0; i < 4; i++)
        d->currscan.comps[i] = 0xFF;
    shdr_comp *el;
    uint8_t comp_pos = 0;
    DL_FOREACH(s->comp, el)
    {
        d->currscan.comps[comp_pos++] = el->cs - 1;
        d->currscan.tbl[el->cs - 1][HUFFMAN_DC] = el->td;
        d->currscan.tbl[el->cs - 1][HUFFMAN_AC] = el->ta;  
    }
    uint8_t fail = 0;
    uint32_t mcuc = 0;
    DL_FOREACH(s->comp, el)
    {
        if (!mcuc)
            mcuc = d->comps[el->cs - 1].mcuc;
        else if (mcuc != d->comps[el->cs - 1].mcuc)
        {
            fail = 1;
            break;
        }
    }
    if (fail)
    {
        free_shdr(s);
        return 0;
    }
    uint16_t ricount = 1;
    d->currscan.ril = mcuc;
    if (d->ril)
    {
        ricount = ceil((double) mcuc / d->ril);
        d->currscan.ril = d->ril;
    }

    d->currscan.ricount = ricount;
    d->currscan.mcur = mcuc;
    d->currscan.mcucurr = 0;
#ifdef DECODER_LOG
    printf("[Decoder]: Restart interval length: %d => Restart intervals: %d\n", d->ril, ricount);
#endif
    fflush(stdout);
    for (uint16_t i = 0; i < ricount; i++)
    {
        if (!decode_restart_interval(d))
            return 0;
    }
    return 1; 
}

// this function decodes an image
static uint8_t decode_frame(decoder_s *d)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding frame\n");
#endif
    uint16_t len;
    fhdr *f = parse_fhdr(d->p, &len);
    if (f->p != 8 || f->nf > 3)
    {
        free_fhdr(f);
        return 0;
    } 
    if (!f->y)
    {
        printf("[Decoder]: This decoder does not support this file\n");
        free_fhdr(f);
        return 0;
    }
    uint8_t fail = 0;
    d->y = f->y;
    d->x = f->x;
    d->nf = f->nf;
    fhdr_comp *el;
    DL_FOREACH(f->comp, el)
    {
        if (el->c > 3)
        {
            fail = 1; 
            break;
        }
        d->comps[el->c - 1].h = el->h;
        d->comps[el->c - 1].v = el->v;
        d->comps[el->c - 1].tq = el->tq;
        if (el->h > d->hmax)
            d->hmax = el->h;
        if (el->v > d->vmax)
            d->vmax = el->v;
    }
    if (fail)
    {
        free_fhdr(f);
        return 0;
    }
    DL_FOREACH(f->comp, el)
    {
        compute_component_params(d, &d->comps[el->c - 1]);
        #ifdef DUMP_FRAME_HEADER
            printf("component %u\n", el->c);
            printf("h %u v %u tq %u\n", el->h, el->v, el->tq);
            component *comp = &d->comps[el->c - 1];
            printf("y %u x %u py %u px %u\n", comp->y, comp->x, comp->py, comp->px);
            printf("blcy %u blcx %u mcuc %u\n", comp->blcy, comp->blcx, comp->mcuc);
        #endif
    }
    free_fhdr(f);
    d->p += len;
    uint16_t marker = get_marker(d);
    while (marker != MRK(EOI))
    {
        while (marker != MRK(SOS))
        {
            if(!interpret_markers(d, marker))
                return 0;
            marker = get_marker(d); 
        }
        if(!decode_scan(d))
            return 0;
#ifdef DECODER_LOG
        printf("[Decoder]: Scan decoded\n");
#endif
        marker = get_marker(d);
    }
    return 1; 
}

// build grayscale image R = G = B = Y
static void build_grayscale(decoder_s *d)
{
    image *img = init_image(d->comps[0].blcy * 8, d->comps[0].blcx * 8);    
    uint32_t bi = d->comps[0].blcy;
    uint32_t bj = d->comps[0].blcx;
    uint32_t bndx = 0;
    for (uint32_t ii = 0; ii < bi; ii++)
    {
        for (uint32_t jj = 0; jj < bj; jj++)
        {
            for (uint8_t i = 0; i < 8; i++)
            {
                for (uint8_t j = 0; j < 8; j++)
                {
                    float p = d->comps[0].m->mat[bndx][i][j];
                    if (p > 255)
                        p = 255;
                    else if (p < 0)
                        p = 0;
                    img->r[ii * 8 + i][jj * 8 + j] = p;
                    img->g[ii * 8 + i][jj * 8 + j] = p;
                    img->b[ii * 8 + i][jj * 8 + j] = p;
                }
            }
            bndx++;
        }
    }
    d->img = img;
}

// extract value corresponding to iimg and jimg image position 
// from a potentially subsampled component
static jpegf get_subsampled_component_value(uint32_t rows, uint32_t cols, 
                                            uint32_t iimg, uint32_t jimg, component *comp)
{
    uint32_t ifactor = rows / (comp->blcy * 8);
    uint32_t jfactor = cols / (comp->blcx  * 8);
    uint32_t icomp = iimg / ifactor;
    uint32_t jcomp = jimg / jfactor;
    uint32_t iblock = icomp / 8;
    uint32_t jblock = jcomp / 8;
    uint32_t iblockpos = icomp % 8;
    uint32_t jblockpos = jcomp % 8;
    uint32_t bndx = iblock * comp->blcx + jblock;
    assert(bndx < comp->blc);
    return comp->m->mat[bndx][iblockpos][jblockpos];
}

// build a color image i.e. convert from YCbCr to RGB
static void build_color(decoder_s *d)
{
    uint32_t rows = d->comps[0].blcy * 8;
    uint32_t cols = d->comps[0].blcx * 8;
    image *img = init_image(rows, cols);    
    uint32_t bi = d->comps[0].blcy;
    uint32_t bj = d->comps[0].blcx;
    uint32_t bndx = 0;
    jpegf y, cb, cr;
    uint8_t r, g, b;
    
    for (uint32_t ii = 0; ii < bi; ii++)
    {
        for (uint32_t jj = 0; jj < bj; jj++)
        {
            // TODO: this can be optimized
            for (uint8_t i = 0; i < 8; i++)
            {
                for (uint8_t j = 0; j < 8; j++)
                {
                    uint32_t iimg = ii * 8 + i;
                    uint32_t jimg = jj * 8 + j;
                    jpegf y = d->comps[0].m->mat[bndx][i][j];
                    jpegf cb = get_subsampled_component_value(rows, cols, iimg, jimg, &d->comps[1]); 
                    jpegf cr = get_subsampled_component_value(rows, cols, iimg, jimg, &d->comps[2]);
                    convert_to_rgb(y, cb, cr, &r, &g, &b);
                    img->r[ii * 8 + i][jj * 8 + j] = r;
                    img->g[ii * 8 + i][jj * 8 + j] = g;
                    img->b[ii * 8 + i][jj * 8 + j] = b;
                }
            } 
            bndx++;
        }
    }
    d->img = img;
}

void *process_component(void *p)
{
    decoder_tp *arg = (decoder_tp *) p;
    decoder_s *d = arg->d;
    uint8_t comp = arg->comp;
    pthread_mutex_t *m = arg->m;
    pthread_cond_t *c = arg->c;
    volatile uint8_t *rdy = arg->rdy;
    
    component *currc = &d->comps[comp];
    int16_t *q = d->tables->quantization[currc->tq];
    uint32_t blc = currc->blc;

    if (blc)
    {
        // allocate the matrices 
        currc->m = alloc_matrices(blc, 8, 8);
        // dequantize, unzigzag, idct and level shift
        for (uint32_t b = 0; b < blc; b++)
        {
            int16_t *currbl = currc->blocks[b];
            jpegf **currmat = currc->m->mat[b];
            dequantize(currbl, q);
            unzigzag(currbl, currmat);
            idct(currmat);
            undo_level_shift(currmat);
        }
    }
    
    pthread_mutex_lock(m);
    *rdy += 1;
    pthread_cond_signal(c); 
    pthread_mutex_unlock(m);
    return NULL;
}

static void obtain_YCbCr(decoder_s *d)
{
    // if only one component no multithreading
    if (d->nf == 1)
    {
        int16_t *q = d->tables->quantization[d->comps[0].tq];
        if (d->comps[0].blc)
        {
            // allocate the matrices 
            d->comps[0].m = alloc_matrices(d->comps[0].blc, 8, 8);
            // dequantize, unzigzag, idct and level shift
            for (uint32_t b = 0; b < d->comps[0].blc; b++)
            {
                int16_t *currbl = d->comps[0].blocks[b];
                jpegf **currmat = d->comps[0].m->mat[b];
                dequantize(currbl, q);
                unzigzag(currbl, currmat);
                idct(currmat);
                undo_level_shift(currmat);
            }
        }
        return;
    }
    // if 3 components do the multithreading, each thread decodes a component
    pthread_mutex_t *m = mutex_create();
    pthread_cond_t *c = cond_create();
    pthread_t *t[3];
    decoder_tp params[3];
    uint8_t rdy = 0;

    for (uint8_t i = 0; i < 3; i++)
    {
        params[i].d = d;
        params[i].comp = i;
        params[i].m = m;
        params[i].c = c;
        params[i].rdy = &rdy;
        t[i] = thread_create(process_component, &params[i]);
    }
    
    pthread_mutex_lock(m);
    while (rdy < 3)
        pthread_cond_wait(c, m);
    pthread_mutex_unlock(m);
    
    for (uint8_t i = 0; i < 3; i++)
        thread_free(t[i]);
    mutex_free(m);
    cond_free(c);
}

static uint8_t rebuild_image(decoder_s *d)
{    
    // obtain the pixel values from the decoded blocks
    obtain_YCbCr(d);
    
    // we won't get rid of padding for now because it is not essential
    // we have YCbCr blocks for all components
    // we have to convert to RGB
    // because are only decoding JFIF according to the standard we have:
    // component 0 => Y
    // component 1 => Cb
    // component 2 => Cr 
    
    // we have grayscale, just copy luma to r g b
    if (d->nf == 1)
        build_grayscale(d);
    else
        build_color(d);
    return 1;
}

// this is the main image decoder function
uint8_t decode_image(decoder *dec)
{
#ifdef DECODER_LOG
    printf("[Decoder]: Decoding image\n");
#endif
    decoder_s *d = (decoder_s *) dec;
    uint8_t *p = d->file->p;
    uint16_t marker = get_marker(d);
    if (marker != MRK(SOI))
        return 0;
    marker = get_marker(d); 
    while (marker != MRK(BDCT_SOF) && marker != MRK(ESDCT_SOF) && marker != MRK(EPDCT_SOF))
    {
        if(!interpret_markers(d, marker))
            return 0;
        marker = get_marker(d);
    }
    uint8_t stat = 0;
    switch (marker)
    {
        case MRK(BDCT_SOF):
        case MRK(ESDCT_SOF):
            d->mode = MODE_BDCT;
            stat = decode_frame(d);
            break;
        default:
            stat = 0;
            break;
    }
    if (stat)
        stat = rebuild_image(d);
    return stat;
}

image *get_image(decoder *dec)
{
    decoder_s *d = (decoder_s *) dec;
    d->gotimg = 1;
    return d->img;
}

// free decoder resources
void free_decoder(decoder *dec)
{
    decoder_s *d = (decoder_s *) dec;
    close_file(d->file);
    if (d->tables)
        free_tables(d->tables);        
    // free matrices
    for (uint8_t i = 0; i < 3; i++)
    {
        if (d->comps[i].blc)
        {
            if (d->comps[i].m)
                free_matrices(d->comps[i].m);
            if (d->comps[i].blocks)
            {
                free(d->comps[i].blocks);
                for (uint32_t b = 0; b < d->comps[i].blc; b++)
                    free(d->comps[i].blocks[b]); 
            }
        }
    }
    if (!d->gotimg && d->img)
        free_image(d->img);
    free(d);
}
