#include <ImageUtil.hpp>

void ImageUtil::resizeImg(const cv::Mat& src, cv::Mat& dst,
                                    int maxSize, bool interpolate) {
  double ratio = 1;
  double w = src.cols;
  double h = src.rows;
  if (w > h)
    ratio = w / (double)maxSize;
  else
    ratio = h / (double)maxSize;
  int nw = (int)(w / ratio);
  int nh = (int)(h / ratio);
  cv::Size sz(nw, nh);
  if (interpolate)
    cv::resize(src, dst, sz);
  else
    resize(src, dst, sz, 0, 0, cv::INTER_NEAREST);
}

void ImageUtil::macOsWaitKey() {
  int i = cv::waitKey();
  cv::destroyAllWindows();
  cv::waitKey(1);
}