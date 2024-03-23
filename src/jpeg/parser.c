#include <jpeg/parser.h>

dht parse_dht(uint8_t *p, uint16_t *l)
{
    uint16_t lh;
    TO_UINT_ADVANCE(uint16_t, lh, p); 
    lh = ENDIAN_SWAP(lh);
    *l = lh;
    lh -= sizeof(uint16_t);
    uint8_t *endp = p + lh;
    dhtt *tables = NULL;
    while (p < endp)
    {
        uint16_t valcount = 0;
        dhtt *t;
        safeMalloc(t, sizeof(dhtt));
        t->tc = *p >> 4;
        t->th = *p & 0xF;
        p++;
        for (uint8_t i = 0; i < 16; valcount += *p, t->l[i++] = *(p++));
        t->lv = valcount;
        safeMalloc(t->v, valcount * sizeof(uint8_t));
        for (uint16_t i = 0; i < valcount; t->v[i++] = *(p++));
        DL_APPEND(tables, t);
    }
    return tables;
}

void free_dht(dht t)
{
    dht el, tmp;
    DL_FOREACH_SAFE(t, el, tmp)
    {
        DL_DELETE(t, el);
        free(el->v);
        free(el);
    }
}

static INLINE void endian_swap_dqt(uint16_t *qt)
{
    for (uint8_t i = 0; i < 64; i++)
        qt[i] = ENDIAN_SWAP(qt[i]);
}

dqt *parse_dqt(uint8_t *p, uint16_t *l)
{
    uint16_t lq;
    TO_UINT_ADVANCE(uint16_t, lq, p);
    lq = ENDIAN_SWAP(lq);
    *l = lq;
    lq -= sizeof(uint16_t);
    dqt *ret = NULL;
    uint8_t *endp = p + lq;
    while (p < endp)
    {
        dqt *t;
        safeMalloc(t, sizeof(dqt));
        t->pq = *p >> 4;
        t->tq = *p & 0xF;
        p++;
        uint16_t qsize = 64 * (sizeof(uint8_t) + t->pq * sizeof(uint8_t));
        safeMalloc(t->q, qsize); 
        memcpy(t->q, p, qsize);
        if (t->pq)
            endian_swap_dqt((uint16_t *) t->q);
        DL_APPEND(ret, t);
        p += qsize;

#ifdef DUMP_QUANTIZATION
        printf("----- QUANTIZATION TABLE -----\n");
        printf("precision %u destination %u\n", t->pq, t->tq);
        if (t->pq)
        {
            uint16_t *table = (uint16_t *) t->q;
            for (uint8_t i = 0; i < 64; i++)
                printf("%u ", table[i]);
            printf("\n");
        }
        else
        {
            uint8_t *table = (uint8_t *) t->q;
            for (uint8_t i = 0; i < 64; i++)
                printf("%u ", table[i]);
            printf("\n");
        }
#endif
    }
   
    return ret;
}

void free_dqt(dqt *q)
{
    dqt *el, *tmp;
    DL_FOREACH_SAFE(q, el, tmp)
    {
        DL_DELETE(q, el);
        free(el->q);
        free(el);
    }
}

static uint16_t __parse_simple(uint8_t *p, uint16_t *l)
{
    uint16_t len;
    TO_UINT_ADVANCE(uint16_t, len, p);
    *l = ENDIAN_SWAP(len);
    uint16_t ret;
    TO_UINT_ADVANCE(uint16_t, ret, p);
    ret = ENDIAN_SWAP(ret);
    return ret;
}

uint16_t parse_dnl(uint8_t *p, uint16_t *l)
{
    return __parse_simple(p, l); 
}

uint16_t parse_dri(uint8_t *p, uint16_t *l)
{
    return __parse_simple(p, l);
}

fhdr *parse_fhdr(uint8_t *p, uint16_t *l)
{
    uint16_t lf;
    TO_UINT_ADVANCE(uint16_t, lf, p);
    lf = ENDIAN_SWAP(lf);
    *l = lf;
    lf -= sizeof(uint16_t);
    uint8_t *endp = p + lf;
    fhdr *ret;
    safeMalloc(ret, sizeof(fhdr));
    TO_UINT_ADVANCE(uint8_t, ret->p, p);
    TO_UINT_ADVANCE(uint16_t, ret->y, p);
    ret->y = ENDIAN_SWAP(ret->y);
    TO_UINT_ADVANCE(uint16_t, ret->x, p);
    ret->x = ENDIAN_SWAP(ret->x);
    TO_UINT_ADVANCE(uint8_t, ret->nf, p);
    ret->comp = NULL;
    for (uint8_t i = 0; i < ret->nf; i++)
    {
        fhdr_comp *c;
        safeMalloc(c, sizeof(fhdr_comp));
        TO_UINT_ADVANCE(uint8_t, c->c, p);
        c->h = *p >> 4;
        c->v = *p & 0xF;
        p++;
        TO_UINT_ADVANCE(uint8_t, c->tq, p);
        DL_APPEND(ret->comp, c);
    }

#ifdef DUMP_FRAME_HEADER
    printf("----- FRAME HEADER -----\n");
    printf("p %u y %u x %u nf %u\n", ret->p, ret->y, ret->x, ret->nf);
    fhdr_comp *el;
    printf("components:\n");
    DL_FOREACH(ret->comp, el)
        printf("h %u v %u tq %u\n", el->h, el->v, el->tq);
#endif

    return ret;
}

void free_fhdr(fhdr *f)
{
    fhdr_comp *el, *tmp;
    DL_FOREACH_SAFE(f->comp, el, tmp)
    {
        DL_DELETE(f->comp, el);
        free(el);
    }
    free(f);
}

shdr *parse_shdr(uint8_t *p, uint16_t *l)
{
    uint16_t ls;
    TO_UINT_ADVANCE(uint16_t, ls, p);
    ls = ENDIAN_SWAP(ls);
    *l = ls;
    ls -= sizeof(uint16_t);
    shdr *ret;
    safeMalloc(ret, sizeof(shdr));
    TO_UINT_ADVANCE(uint8_t, ret->ns, p);
    ret->comp = NULL;
    for (uint8_t i = 0; i < ret->ns; i++)
    {
        shdr_comp *c;
        safeMalloc(c, sizeof(shdr_comp));
        TO_UINT_ADVANCE(uint8_t, c->cs, p);
        c->td = *p >> 4;
        c->ta = *p & 0xF;
        p++;
        DL_APPEND(ret->comp, c);
    }
    TO_UINT_ADVANCE(uint8_t, ret->ss, p);
    TO_UINT_ADVANCE(uint8_t, ret->se, p);
    ret->ah = *p >> 4;
    ret->al = *p & 0xF;

#ifdef DUMP_SCAN_HEADER
    printf("----- SCAN HEADER -----\n");
    printf("ns: %u\n", ret->ns);
    shdr_comp *el;
    DL_FOREACH(ret->comp, el)
        printf("cs %u td %u ta %u\n", el->cs, el->td, el->ta);
    printf("ss %u se %u ah %u al %u\n", ret->ss, ret->se, ret->ah, ret->al);
#endif
    
    return ret;
}

void free_shdr(shdr *s)
{
    shdr_comp *el, *tmp;
    DL_FOREACH_SAFE(s->comp, el, tmp)
    {
        DL_DELETE(s->comp, el);
        free(el);
    }
    free(s);
}

uint16_t get_marker_seg_len(uint8_t *p)
{
    uint16_t ret;
    TO_UINT_ADVANCE(uint16_t, ret, p);
    return ENDIAN_SWAP(ret);
}
