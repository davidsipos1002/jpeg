#include <parser.h>

dht parse_dht(uint8_t *p)
{
    uint16_t lh;
    TO_UINT_ADVANCE(uint16_t, lh, p); 
    lh = ENDIAN_SWAP(lh) - sizeof(uint16_t);
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

dqt *parse_dqt(uint8_t *p)
{
    uint16_t lq;
    TO_UINT_ADVANCE(uint16_t, lq, p);
    lq = ENDIAN_SWAP(lq) - sizeof(uint16_t);
   
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
