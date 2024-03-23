#include <iostream>

#include <file.h>
#include <markers.h>
#include <parser.h>
#include <huffman.h>

int main()
{
    mmapfile *f = map_file("img/bdct_simple.jpg");
    if (!f)
    {
        std::cout << "fail";
        return 0;
    }
    uint8_t *p = (uint8_t *) f->p;
    for (size_t i = 0; i < f->size; i++)
    {
        uint16_t *marker = (uint16_t *) p;
        if (ENDIAN_SWAP(*marker) == MRK(SOI))
            std::cout << "start of image" << std::endl;
        else if (ENDIAN_SWAP(*marker) == MRK(EOI))
            std::cout << "end of image" << std::endl;
        else if (ENDIAN_SWAP(*marker) == MRK(DHT)) 
        {
            std::cout << "define huffman table" << std::endl;
            
            dht tables = parse_dht(p + 2);
            dht curr;
            DL_FOREACH(tables, curr)
            {
                extract_huffman_table(curr);
            }
            free_dht(tables);
        }
        else if (ENDIAN_SWAP(*marker) == MRK(BDCT_SOF))
            std::cout << "baseline sof" << std::endl;
        p++;
    }
    close_file(f);
    return 0;
}