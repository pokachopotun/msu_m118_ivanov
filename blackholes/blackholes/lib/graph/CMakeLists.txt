cmake_minimum_required(VERSION 3.2)

project(graph)

add_library(graph graph.cpp)

set_target_properties(graph PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(graph
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../common"
)

if (MSVC)
    target_compile_options(graph PRIVATE
        /W4
    )
else ()
    target_compile_options(graph PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()
