cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 23)

find_package(OpenCV REQUIRED)

project(JPEG VERSION 0.1 LANGUAGES CXX)

set(includes ${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/nativefiledialog-extended)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/external/util)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/imgjpeg)

add_executable(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${OpenCV_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC ${OpenCV_LIBS})
target_link_libraries(${PROJECT_NAME} PUBLIC nfd)
target_link_libraries(${PROJECT_NAME} PUBLIC util)
target_link_libraries(${PROJECT_NAME} PUBLIC imgjpeg)
