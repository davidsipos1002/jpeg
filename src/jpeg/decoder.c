#include <jpeg/decoder.h>

#include <assert.h>
#include <math.h>

#include <file.h>
#include <jpeg/tables.h>
#include <jpeg/parser.h>

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
    uint8_t **blocks; // decoded blocks of the component
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
} decoder_s;

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
    safeMalloc(c->blocks, block_count * sizeof(uint8_t *));
    for (uint32_t i = 0; i < block_count; i++)
    {
        safeMalloc(c->blocks[i], 64 * sizeof(uint8_t));
        memset(c->blocks[i], 0, 64 * sizeof(uint8_t));
    }
    c->blc = block_count;
    c->mcuc = (c->blcy / c->v) * (c->blcx / c->h);
}

decoder *init_decoder(const char *filename)
{
    decoder_s *ret;
    memset(ret, 0, sizeof(decoder_s));
    safeMalloc(ret, sizeof(decoder_s));
    ret->file = map_file(filename);
    ret->p = ret->file->p;
    ret->tables = init_tables();
    return (decoder *) ret;
}

static uint16_t get_marker(decoder_s *d)
{
    uint16_t marker;
    TO_UINT_ADVANCE(uint16_t, marker, d->p); 
    return ENDIAN_SWAP(marker);
}

static void interpret_markers(decoder_s *d, uint16_t marker)
{
    uint16_t len;
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
        default:
        {
            len = get_marker_seg_len(d->p);
            break;
        }
    }
    d->p += len;
}

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

static uint8_t nextbit(decoder_s *d)
{
    static uint8_t b;
    if (!d->currscan.cnt)
    {
        TO_UINT_ADVANCE(uint8_t, b, d->p);
        d->currscan.cnt = 8;
        if (b == 0xFF)
        {
            uint8_t b2 = TO_UINT_ADVANCE(uint8_t, b2, d->p); 
            assert(b2 == 0);
        }
    }
    uint8_t ret = b >> 7;
    d->currscan.cnt--;
    b <<= 1;
    return ret;
}

static int16_t receive(decoder_s *d, uint8_t s)
{
    int16_t v = 0;
    for (uint8_t i = 0;i < s; i++)
        v = (v << 1) + nextbit(d);
    return v;
}

uint8_t decode(decoder_s *d, uint8_t coeff)
{
    uint8_t i = 0;
    uint8_t code = nextbit(d);
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
    return hft->huffval[j];
}

static uint8_t decode_data_unit(decoder_s *d, uint8_t *block)
{
    // decode DC DIFF
    uint8_t t = decode(d, HUFFMAN_DC);
    int16_t diff = receive(d, t);
    diff = extend(diff, t);
    // obtain DC coefficient
    block[0] = diff + d->currscan.pred[d->currscan.currc];
    d->currscan.pred[d->currscan.currc] = block[0];
    
    // decode AC coefficients
    uint8_t k = 1;
    while (k < 64)
    {
        uint8_t rs = decode(d, HUFFMAN_AC);
        uint8_t s = rs & 0xF;
        uint8_t r = rs >> 4;
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
        block[k] = extend(block[k], s);
        k++;
    }
    return 1;
}

static uint8_t decode_mcu(decoder_s *d)
{
    uint8_t cpos = 0;
    uint8_t c;
    while ((c = d->currscan.comps[cpos]) != 0xFF)
    {
        component *comp = &d->comps[c]; 
        uint16_t mcu_rows = comp->blcy / comp->v;
        uint16_t mcu_cols = comp->blcx / comp->h;
        uint16_t mcuy = d->currscan.mcucurr / mcu_cols;
        d->currscan.currc = c;
        if (mcuy >= mcu_rows)
            return 0;
        uint16_t mcux = d->currscan.mcucurr % mcu_cols; 
        for (uint8_t i = 0; i < comp->v; i++)
        {
            for (uint8_t j = 0; j < comp->h; j++)
            {
                uint16_t bi = mcuy * comp->v + i;
                uint16_t bj = mcux * comp->h + j;
                uint8_t *block = comp->blocks[bi * comp->blcx + bj];
                decode_data_unit(d, block);
            }
        }
        cpos++; 
    }
    d->currscan.mcucurr++;
    return 1;
}

static uint8_t is_marker(uint16_t marker)
{
    return (marker >> 8) == 0xFF && (marker & 0xFF) != 0;
}

static uint8_t decode_restart_interval(decoder_s *d)
{
    // reset decoder
    memset(d->currscan.pred, 0, 3 * sizeof(uint8_t));
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

static uint8_t decode_scan(decoder_s *d)
{
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
#ifdef DUMP_SCAN_HEADER
    printf("ril: %d => restart intervals: %d\n", d->ril, ricount);
#endif
    fflush(stdout);
    for (uint16_t i = 0; i < ricount; i++)
    {
        if (!decode_restart_interval(d))
            return 0;
    }
    return 1; 
}

static uint8_t decode_frame(decoder_s *d)
{
    uint16_t len;
    fhdr *f = parse_fhdr(d->p, &len);
    if (f->p != 8 || f->nf > 3)
    {
        free_fhdr(f);
        return 0;
    } 
    if (!f->y)
    {
        printf("this decoder does not support this file\n");
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
        #ifdef DUMP_COMPONENT_PARAMS
            printf("component %u\n", el->c);
            printf("h %u v %u tq %u\n", el->h, el->v, el->tq);
            component *comp = &d->comps[el->c - 1];
            printf("y %u x %u py %u px %u\n", comp->y, comp->x, comp->py, comp->px);
            printf("blcy %U blcx %u mcuc %u\n", comp->blcy, comp->blcx, comp->mcuc);
        #endif
    }
    free_fhdr(f);
    d->p += len;
    uint16_t marker = get_marker(d);
    while (marker != MRK(EOI))
    {
        while (marker != MRK(SOS))
        {
            interpret_markers(d, marker);
            marker = get_marker(d); 
        }
        if(!decode_scan(d))
            return 0;
        marker = get_marker(d);
    }
    return 1; 
}

uint8_t decode_image(decoder *dec)
{
    decoder_s *d = (decoder_s *) dec;
    uint8_t *p = d->file->p;
    uint16_t marker = get_marker(d);
    if (marker != MRK(SOI))
        return 0;
    marker = get_marker(d); 
    while (marker != MRK(BDCT_SOF) && marker != MRK(ESDCT_SOF) && marker != MRK(EPDCT_SOF))
    {
        interpret_markers(d, marker);
        marker = get_marker(d);
    }
    switch (marker)
    {
        case MRK(BDCT_SOF):
            d->mode = MODE_BDCT;
            return decode_frame(d); 
    }
    return 0;
}

void free_decoder(decoder *dec)
{
    decoder_s *d = (decoder_s *) dec;
    close_file(d->file);
    free_tables(d->tables);        
}
