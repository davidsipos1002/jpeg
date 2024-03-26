cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

project(util VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_GENERATOR "ninja")

find_package(OpenCV REQUIRED)

set(includes ${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.c)

add_library(${PROJECT_NAME} STATIC ${sources})

target_compile_definitions(${PROJECT_NAME} PUBLIC ASSETS_DIR="${CMAKE_SOURCE_DIR}/external/assets")
target_compile_definitions(${PROJECT_NAME} PUBLIC IMAGE_DIR="${CMAKE_SOURCE_DIR}/external/assets/Images")
target_compile_definitions(${PROJECT_NAME} PUBLIC VIDEO_DIR="${CMAKE_SOURCE_DIR}/external/assets/Videos")
target_include_directories(${PROJECT_NAME} PUBLIC ${includes} ${OpenCV_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC nfd ${OpenCV_LIBS})
