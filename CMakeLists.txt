cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 23)

project(JPEG VERSION 0.1 LANGUAGES C CXX)

set(includes ${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.c ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes})
target_link_libraries(${PROJECT_NAME} PUBLIC pthread)
# target_compile_definitions(${PROJECT_NAME} PRIVATE DECODER_LOG DUMP_HUFFMAN DUMP_QUANTIZATION DUMP_FRAME_HEADER DUMP_SCAN_HEADER)
