#include <mem.h>

#include <string.h>

arrmat *alloc_matrices(uint32_t count, uint32_t row, uint32_t col)
{
    uint32_t matsize = count * row * col * sizeof(float);
    float *store;
    safeMalloc(store, matsize);
    memset(store, 0, matsize);
    float ***mats;
    safeMalloc(mats, count * sizeof(float **));
    for (uint32_t i = 0; i < count; i++)
        safeMalloc(mats[i], row * sizeof(float *));
    for (uint32_t i = 0; i < count; i++)
        for (uint32_t j = 0; j < row; j++)
            mats[i][j] = &store[i * row * col + j * row]; 
    arrmat *ret;
    safeMalloc(ret, sizeof(arrmat));
    ret->count = count;
    ret->row = row;
    ret->col = col;
    ret->mat = mats;
    return ret;
}

void free_matrices(arrmat *p)
{
    float *store = p->mat[0][0];
    free(store);
    for (uint32_t i = 0; i < p->count; i++)
        free(p->mat[i]);
    free(p->mat);
    free(p);
}
