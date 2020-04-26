cmake_minimum_required(VERSION 3.2)

project(algorithm)

add_library(algorithm algorithm.cpp)

set_target_properties(algorithm PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(algorithm
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../common"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../graph"
)

if (MSVC)
    target_compile_options(algorithm PRIVATE
        /W4
    )
else ()
    target_compile_options(algorithm PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()
