cmake_minimum_required(VERSION 3.20)

project(lyapunov LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

if(NOT (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
    find_package(OpenMP REQUIRED)
endif()

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")

if(WIN32 AND NOT MSVC)
    # Statically link GCC's standard libraries, but allow shared runtimes such as OpenMP.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")

    # Prefer static variations (.a) over dynamic (.dll.a) when finding third-party libraries
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(BUILD_SHARED_LIBS OFF)
endif()

add_executable(lyapunov lyapunov.cc
        fractal_type.cc
        lissajous.cc
        string_case.hh)

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(lyapunov PRIVATE -fopenmp)
    target_link_options(lyapunov PRIVATE -fopenmp)
else()
    target_link_libraries(lyapunov PRIVATE OpenMP::OpenMP_CXX)
endif()

add_executable(gen-ab-sequences gen-ab-sequences.cc)

add_executable(gen_zoom gen_zoom.cc)

add_executable(gen_poses gen_poses.cc
        fractal_type.cc)

add_executable(lissa_overview lissajous.cc lissa_overview.cc)

install(TARGETS lyapunov gen_zoom gen-ab-sequences gen_poses lissa_overview
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
