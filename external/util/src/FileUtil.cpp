#include <FileUtil.hpp>

std::string FileUtil::getFileOrDir(const bool isFile) {
  NFD::Guard nfdGuard;
  NFD::UniquePath outPath;
  nfdresult_t result;
  if (isFile) {
    result = NFD::OpenDialog(outPath);
  } else {
    result = NFD::PickFolder(outPath);
  }
  if (result == NFD_OKAY) {
    return outPath.get();
  }
  return "";
}

std::string FileUtil::getSingleFileAbsPath() {
  return getFileOrDir();
}

std::string FileUtil::getDirectoryAbsPath() {
  return getFileOrDir(false);
}

std::vector<std::string> FileUtil::getAllFilesInDirectory() {
  std::vector<std::string> filenames;
  const std::string dir_abs_path = getDirectoryAbsPath();
  if (!dir_abs_path.empty()) {
    const stdfs::path path_to_traverse(dir_abs_path);
    if (stdfs::exists(path_to_traverse) && stdfs::is_directory(path_to_traverse)) {
      for (const auto &entry : stdfs::directory_iterator(path_to_traverse)) {
        auto abs_file_path = entry.path().string();
        filenames.push_back(abs_file_path);
      }
    }
  }
  return filenames;
}