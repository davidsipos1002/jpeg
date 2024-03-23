#include <jpeg/tables.h>

tables *init_tables()
{
    tables *t;
    safeMalloc(t, sizeof(tables));
    memset(t, 0, sizeof(tables));
    return t;
}

void install_huffman_table(tables *t, hfft *h)
{
    hfft **ts = (h->tclass == HUFFMAN_DC) ? t->huffman_dc : t->huffman_ac;
    if (ts[h->destination])
        free_huffman_table(ts[h->destination]);
    ts[h->destination] = h;
}

void install_quantization_table(tables *t, dqt *q)
{
    if (q->pq == QUANTIZATION_16)
    {
        printf("only 8-bit images are supported\n");
        exit(1);
    }
    uint8_t **ts = t->quantization;
    if (ts[q->tq])
        free(ts[q->tq]);
    safeMalloc(ts[q->tq], 64 * sizeof(uint8_t));
    memcpy(ts[q->tq], q->q, 64 * sizeof(uint8_t));
}

void free_tables(tables *t)
{
    for (uint8_t i = 0; i < 4; i++)
    {
        if (t->huffman_dc[i])
            free_huffman_table(t->huffman_dc[i]);
        if (t->huffman_ac[i])
            free_huffman_table(t->huffman_ac[i]);
        if (t->quantization[i])
            free(t->quantization[i]);
    }
}
