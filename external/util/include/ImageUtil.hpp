#pragma once

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

class ImageUtil {
 public:
  static void resizeImg(const cv::Mat& src, cv::Mat& dst, int maxSize,
                        bool interpolate);
  static void macOsWaitKey();
};
