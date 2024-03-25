#include <jpeg/image.h>

#include <mem.h>

image *init_image(uint16_t rows, uint16_t cols)
{
    image *ret;
    safeMalloc(ret, sizeof(image));
    ret->y = rows;
    ret->x = cols;
    ret->r = alloc_mat(rows, cols);
    ret->g = alloc_mat(rows, cols);
    ret->b = alloc_mat(rows, cols);
    return ret;
}

void free_image(image *img)
{
    free_mat(img->r, img->y, img->x);
    free_mat(img->g, img->y, img->x);
    free_mat(img->b, img->y, img->x);
    free(img);
}
