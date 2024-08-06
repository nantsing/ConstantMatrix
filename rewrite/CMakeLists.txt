cmake_minimum_required(VERSION 3.12)
project(MatrixMul)

set(CMAKE_CXX_STANDARD 20)

# Add your source files here
set(SOURCES
    main.cpp
)

# Add your header files here
set(HEADERS
    csd.hpp
    adder.hpp
    circuit.hpp
)

add_executable(MatrixMul ${SOURCES} ${HEADERS})