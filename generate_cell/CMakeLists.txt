cmake_minimum_required(VERSION 3.16)
project(generate_cell LANGUAGES CXX)

# set compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -march=native")

# project source files
file(GLOB SRCFILES
        src/*.cpp
        src/Mesh/*.cpp
        src/Utils/*.cpp
)

# remove main.cpp from ${SRCFILES}
FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "src/main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "src/main.cpp")
ENDFOREACH(item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC
        src
        src/Mesh
        src/Utils
)

add_executable(${PROJECT_NAME}_bin src/main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)

# code only support 3D
message(STATUS "Building for three dimensions (3D)")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC DIM_=3)

# add -pthread to compilation and linking
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${CMAKE_THREAD_LIBS_INIT})