#pragma once

#include <filesystem>
#include <iostream>
#include <numeric>
#include <vector>

#include <nfd.hpp>

namespace stdfs = std::filesystem;

class FileUtil {
 private:
  static std::string getFileOrDir(bool isFile = true);

 public:
  static std::string getSingleFileAbsPath();
  static std::string getDirectoryAbsPath();
  static std::vector<std::string> getAllFilesInDirectory();
};
