cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 23)

project(jpeg VERSION 0.1 LANGUAGES C CXX)

set(includes ${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB_RECURSE sources ${CMAKE_CURRENT_LIST_DIR}/src/*.c ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes})

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_LIST_DIR})