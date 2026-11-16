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

# The library with the fractal computation, usable by the CLI and by GUI programs.
add_library(lyapunovlib STATIC
        fractal_params.cc
        fractal_type.cc
        grid.cc
        file_writer.cc
        render.cc
        draw-rgb.cc
        fonts/fonts.cc
        fractal_params.hh
        fractal_type.hh
        grid.hh
        file_writer.hh
        render.hh
        render.hxx
        draw.hh
        draw.hxx
        rgb.hh
        string_case.hh)

add_executable(lyapunov lyapunov.cc)
target_link_libraries(lyapunov PRIVATE lyapunovlib)

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(lyapunovlib PRIVATE -fopenmp)
    target_link_options(lyapunov PRIVATE -fopenmp)
else()
    target_link_libraries(lyapunovlib PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(gen-ab-sequences gen-ab-sequences.cc)

add_executable(gen_zoom gen_zoom.cc)

add_executable(gen_poses gen_poses.cc)
target_link_libraries(gen_poses PRIVATE lyapunovlib)

add_executable(generate_font generate_font.cc)

add_executable(lissa_overview lissajous.cc lissa_overview.cc)

install(TARGETS lyapunov gen_zoom gen-ab-sequences gen_poses generate_font lissa_overview
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# The library for external consumers (e.g. GUI programs): static archive
# plus the public headers (the *.hxx template definitions are internal).
install(TARGETS lyapunovlib ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES fractal_params.hh fractal_type.hh grid.hh file_writer.hh
        render.hh draw.hh rgb.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lyapunov)
