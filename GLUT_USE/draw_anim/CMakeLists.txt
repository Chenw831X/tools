cmake_minimum_required(VERSION 3.16)
project(draw_anim LANGUAGES CXX)

# include header file of Eigen
include_directories(/usr/local/include)

add_executable(${PROJECT_NAME} main.cpp)

find_package(GLUT REQUIRED) # GLUT
target_link_libraries(${PROJECT_NAME} PUBLIC glut)

find_package(OpenGL REQUIRED) # opengl
target_link_libraries(${PROJECT_NAME} PUBLIC OpenGL)

# add -pthread to compilation and linking
#find_package(Threads)
#target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_THREAD_LIBS_INIT})