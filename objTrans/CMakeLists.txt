cmake_minimum_required(VERSION 3.16)
project(objTrans LANGUAGES CXX)

add_executable(${PROJECT_NAME} src/main.cpp)

# add -pthread to compilation and linking
find_package(Threads)
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_THREAD_LIBS_INIT})