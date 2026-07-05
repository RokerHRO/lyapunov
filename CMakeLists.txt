cmake_minimum_required(VERSION 3.10)

project(lyapunov LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O0 -ggdb -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ggdb -fopenmp")

add_executable(lyapunov lyapunov.cc)

install(TARGETS lyapunov
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
