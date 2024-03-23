#pragma once

#include <stdlib.h>
#include <stdint.h>

#include <misc.h>

BEGIN_C_DECLS

typedef struct
{
    size_t size;
    void *p;
} mmapfile;

mmapfile *map_file(const char *filename);
void close_file(mmapfile *file);

END_C_DECLS