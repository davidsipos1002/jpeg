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