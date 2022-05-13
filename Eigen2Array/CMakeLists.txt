cmake_minimum_required(VERSION 3.9)
project(Eigen2Array LANGUAGES CXX)

add_executable(${PROJECT_NAME} src/main.cpp)
include_directories("/usr/local/include/eigen-3.4.0")