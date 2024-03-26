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
    // run the decoder
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

void encodeImage(const std::string &source, const std::string &destination, uint8_t quality)
{
    // read image
    cv::Mat src = cv::imread(source, cv::IMREAD_COLOR);

    // build the image structure for the encoder
    image *img = init_image(src.rows, src.cols);
    for (uint16_t i = 0; i < src.rows; i++)
    {
        for (uint16_t j = 0; j < src.cols; j++)
        {
            const cv::Vec3b &p = src.at<cv::Vec3b>(i, j);
            img->b[i][j] = p[0];
            img->g[i][j] = p[1];
            img->r[i][j] = p[2];
        }
    }
    
    // run the encoder
    uint8_t q = 0;
    switch (quality)
    {
        case 100:
            q = JPEG_QUALITY_100;
            break;
        case 80:
            q = JPEG_QUALITY_80;
            break;
        case 50:
            q = JPEG_QUALITY_50;
            break;
    }
    encoder *e = init_encoder(destination.c_str(), img, q);
    uint8_t stat = encode_image(e);
    free_encoder(e);
    if (!stat)
        std::cout << "Failed to encode image!" << std::endl;
    free_image(img);
}

int main()
{
    char c;
    while (true)
    {
        std::cout << "Select between d (decode), e (encode), x (exit) : ";
        std::cin >> c;
        
        if (c == 'x')
            break;
        std::cout << "Select source file" << std::endl;
        std::string src = FileUtil::getSingleFileAbsPath();
        if (src.empty())
        {
            std::cout << "File not selected!" << std::endl;
            continue;
        }
        std::filesystem::path path(src);
        std::string dst = path.filename().string();
        dst = c == 'e' ? dst + "-encoded.jpg" : dst + ".dump";
        if (c == 'd')
        {
            std::cout << "Do you want to dump the pixel values ? (y/n) ";
            std::cin >> c;
            if (c == 'n')
                dst = "";
            decodeImage(src, dst);
        }
        else
        {
            std::cout << "Enter quality (50, 80, 100): ";
            int q;
            std::cin >> q;
            if (q != 50 && q != 80 && q != 100)
            {
                std::cout << "Invalid quality!" << std::endl;
                continue;
            }
            encodeImage(src, dst, static_cast<uint8_t>(q));
        }
    }
    return 0;
}
