cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

project(imgjpeg VERSION 0.1 LANGUAGES C)

file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.c ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_library(${PROJECT_NAME} STATIC ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)

target_link_libraries(${PROJECT_NAME} PRIVATE pthread)

# target_compile_definitions(${PROJECT_NAME} PRIVATE DECODER_LOG DUMP_HUFFMAN DUMP_QUANTIZATION DUMP_FRAME_HEADER DUMP_SCAN_HEADER)
