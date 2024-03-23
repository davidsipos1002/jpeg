#include <file.h>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

mmapfile *map_file(const char *filename)
{
    int fd = open(filename, O_RDONLY);
    if (fd == -1)
        return NULL;
    size_t sz = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, SEEK_SET);
    void *p = mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fd, 0); 
    if (p == MAP_FAILED)
        return NULL;
    mmapfile *ret;
    safeMalloc(ret, sizeof(mmapfile));
    ret->size = sz;
    ret->p = p;
    return ret; 
}

void close_file(mmapfile *file)
{
    munmap(file->p, file->size);
    free(file);
}