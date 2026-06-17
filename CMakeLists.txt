cmake_minimum_required(VERSION 3.5)

project(lyapunov LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(lyapunov lyapunov.cc)

install(TARGETS lyapunov
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
