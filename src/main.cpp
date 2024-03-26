#include <iostream>
#include <fstream>
#include <filesystem>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/opencv.hpp>

#include <FileUtil.hpp>
#include <ImageUtil.hpp>

#include <jpeg/image.h>
#include <jpeg/decoder.h>
#include <jpeg/encoder.h>

void dumpComponent(std::ofstream &out, uint16_t rows, uint16_t cols, uint8_t **comp)
{
    for (uint16_t i = 0; i < rows; i++)
    {
        for (uint16_t j = 0; j < cols; j++)
            out << (uint32_t) comp[i][j] << " ";
        out << "\n";
    }
}

void decodeImage(const std::string &source, const std::string &destination)
{
    decoder *d = init_decoder(source.c_str());
    if (!decode_image(d))
    {
        free_decoder(d);
        std::cout << "Failed to decode image!" << std::endl;
        return;
    }
    image *img = get_image(d);
    free_decoder(d);

    cv::Mat src(img->y, img->x, CV_8UC3, cv::Scalar(0, 0, 0));
    for (int i = 0; i < img->y; i++)
    {
        for (int j = 0; j < img->x; j++)
        {
            cv::Vec3b &p = src.at<cv::Vec3b>(i, j);
            p[0] = img->b[i][j];
            p[1] = img->g[i][j];
            p[2] = img->r[i][j];
        }
    }

    cv::imshow("Decoded image", src);
    ImageUtil::macOsWaitKey();

    if (!destination.empty())
    {
        std::cout << "Dumping pixel values..." << std::endl;
        std::ofstream out(destination);

        out << "Dimensions: " << img->y << " " << img->x << "\n\n";

        out << "R\n";
        dumpComponent(out, img->y, img->x, img->r); 
        out << "\n";

        out << "G\n";
        dumpComponent(out, img->y, img->x, img->r); 
        out << "\n";
        
        out << "B\n";
        dumpComponent(out, img->y, img->x, img->r);

        out.close();
        
        std::cout << "Done." << std::endl;
    }
    
    free_image(img);
}

void encodeImage()
{
    decoder *d = init_decoder("img/everest_sub.jpg");
    decode_image(d);
    image *img = get_image(d);
    free_decoder(d);
    encoder *e = init_encoder("test.jpg", img, JPEG_QUALITY_100);
    encode_image(e);
    free_encoder(e);
}

int main()
{
    // char c;
    // while (true)
    // {
    //     std::cout << "Select between d (decode), e (encode), x (exit) : ";
    //     std::cin >> c;
        
    //     if (c == 'x')
    //         break;
    //     if (c == 'e')
    //     {
    //         std::cout << "Encoding not supported!" << std::endl;
    //         continue;
    //     }

    //     std::cout << "Select source file for decoding" << std::endl;
    //     std::string src = FileUtil::getSingleFileAbsPath();
    //     if (src.empty())
    //     {
    //         std::cout << "File not selected!" << std::endl;
    //          continue;
    //     }
    //     std::cout << "Do you want to dump the pixel values ? (y/n) ";
    //     std::cin >> c;
    //     std::filesystem::path path(src);
    //     std::string dst = c == 'y' ? path.filename().string() + ".dump" : "";
    //     decodeImage(src, dst);
    // }
    encodeImage();
    return 0;
}
