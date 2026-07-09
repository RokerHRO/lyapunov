cmake_minimum_required(VERSION 3.10)

project(lyapunov LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O0 -ggdb -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ggdb -fopenmp")

if(WIN32 AND NOT MSVC)
    # Force static linking of standard libraries
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")

    # Prefer static variations (.a) over dynamic (.dll.a) when finding third-party libraries
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
endif()

add_executable(lyapunov lyapunov.cc)
add_executable(gen-ab-sequences gen-ab-sequences.cc)
add_executable(gen_zoom gen_zoom.cc)

install(TARGETS lyapunov
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
