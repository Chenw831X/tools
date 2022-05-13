cmake_minimum_required(VERSION 3.16)
project(LinSysSolver LANGUAGES C CXX)

# alternative: CHOLMOD, EIGEN
set(LINSYSSOLVER "CHOLMOD" CACHE STRING "Linear system solver to use (options: CHOLMOD, EIGEN)")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# set compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -march=native -m64")

# project source files
file(GLOB SRCFILES src/LinSysSolver/*.cpp)
add_library(${PROJECT_NAME}_dev ${SRCFILES})
target_include_directories(${PROJECT_NAME}_dev PUBLIC src/LinSysSolver)

add_executable(${PROJECT_NAME}_bin src/main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)

# tbb
find_package(TBB REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC TBB::tbb)
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_TBB)

if(LINSYSSOLVER STREQUAL "CHOLMOD")
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_CHOLMOD)
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC USE_CHOLMOD)
else()
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC LINSYSSOLVER_USE_EIGEN)
endif()

# code only support 3D
message(STATUS "Building for three dimensions (3D)")
target_compile_definitions(${PROJECT_NAME}_dev PUBLIC DIM=3)

# add -pthread to compilation and linking
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${CMAKE_THREAD_LIBS_INIT})

# MKL
set(MKL_INTERFACE "lp64")
#set(MKL_THREADING "gnu_thread")
#set(ENABLE_BLAS95 ON)
#set(ENABLE_LAPACK95 ON)
find_package(MKL CONFIG REQUIRED)
if(MKL_FOUND)
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC MKL::MKL)
    #target_include_directories(${PROJECT_NAME}_dev PUBLIC ${MKL_INCLUDE})
    target_compile_definitions(${PROJECT_NAME}_dev PUBLIC EIGEN_USE_MKL_ALL)
endif()


if(LINSYSSOLVER STREQUAL "CHOLMOD")
    # SuiteSparse
    find_package(SuiteSparse REQUIRED)

    target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_LIBRARIES})
    target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SUITESPARSE_INCLUDE_DIRS})
endif()
