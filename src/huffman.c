#include <huffman.h>

static uint8_t *gen_huffsize(dht t)
{
    uint8_t *sztable;
    safeMalloc(sztable, t->lv * sizeof(uint8_t));
    uint16_t k = 0;
    for (uint8_t i = 0; i < 16; i++)
        for (uint8_t j = 0; j < t->l[i]; j++)
            sztable[k++] = i + 1;
   return sztable;
}

static uint16_t *gen_huffcode(uint8_t *huffsize, uint16_t lastk)
{
    uint16_t *huffcode;
    safeMalloc(huffcode, lastk * sizeof(uint16_t));
    uint16_t code = 0;
    uint16_t k = 0;
    uint8_t currsz = huffsize[0];
    while (k < lastk)
    {
        if (huffsize[k] == currsz)
        {
            huffcode[k] = code;
            code++;
            k++;
        }
        else
        {
            while (huffsize[k] != currsz)
            {
                code <<= 1;
                currsz++;
            }
        }
    }
    return huffcode;
}

static void order_codes(dht t,uint8_t *huffsize, uint16_t *huffcode, hfft *huff)
{
    memset(huff->ehufco, 0, 256 * sizeof(uint16_t));
    memset(huff->ehufsi, 0, 256 * sizeof(uint8_t));
    for (uint16_t i = 0; i < t->lv; i++)
    {
        uint8_t val = t->v[i];
        huff->ehufco[val] = huffcode[i];
        huff->ehufsi[val] = huffsize[i];
    }
}

static void gen_decoder_tables(dht t, hfft *huff, uint16_t *huffcode)
{
    uint16_t j = 0;
    for (uint8_t i = 0; i < 16; i++)
    {
        if(!t->l[i])
        {
            huff->maxcode[i] = 0xFFFF;
            continue;
        }
        huff->valptr[i] = j;
        huff->mincode[i] = huffcode[j];
        j += t->l[i] - 1;
        huff->maxcode[i] = huffcode[j];
        j++;
    }
}

hfft *extract_huffman_table(dht table)
{
    uint8_t *huffsize = gen_huffsize(table);
    uint16_t *huffcode = gen_huffcode(huffsize, table->lv);
    hfft *huff;
    safeMalloc(huff, sizeof(hfft));
    huff->valcount = table->lv;
    safeMalloc(huff->huffval, table->lv * sizeof(uint8_t));
    order_codes(table, huffsize, huffcode, huff); 
    memcpy(huff->huffval, table->v, table->lv * sizeof(uint8_t));
    gen_decoder_tables(table, huff, huffcode);

#ifdef DUMP_HUFFMAN
    printf("category: %u destination: %u\n", table->tc, table->th);
    printf("bits:\n");
    for (uint16_t i = 0;i < 16;i ++)
       printf("%u ", table->l[i]);
    printf("\n");
    printf("huffval:\n");
    for (uint16_t i = 0;i < table->lv; i++)
       printf("%u ", table->v[i]);
    printf("\n");
    printf("huffsize:\n");
    for (uint16_t i = 0; i < table->lv; i++)
        printf("%d ", huffsize[i]);
    printf("\n");
    printf("huffcode:\n");
    for (uint16_t i = 0; i < table->lv; i++)
    {
        print_binary16(huffcode[i]);
        printf("\n");
    }
    printf("encoder table:\n");
    for (uint16_t i = 0; i < 256; i++)
    {
        printf("value %u => size: %u code: ", i, huff->ehufsi[i]);
        print_binary16(huff->ehufco[i]);
        printf("\n");
    }
    
    printf("decoder table:\n");
    for (uint8_t i = 0; i < 16; i++)
        printf("min %u max %u valptr %u\n", huff->mincode[i], huff->maxcode[i], huff->valptr[i]);
#endif

    free(huffsize);
    free(huffcode);
    return NULL;
}