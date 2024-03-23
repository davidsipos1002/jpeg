#include <iostream>

#include <jpeg/decoder.h>

int main()
{
    decoder *dec = init_decoder("img/lena_small.jpg");
    uint32_t ret =  decode_image(dec);
    std::cout << "Result: " << ret << std::endl;
    free_decoder(dec);
    return 0;
}