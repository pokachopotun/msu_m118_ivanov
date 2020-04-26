cmake_minimum_required(VERSION 3.2)

project(io)

add_library(io io.cpp)

set_target_properties(io PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(io
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../common"
)

if (MSVC)
    target_compile_options(io PRIVATE
        /W4
    )
else ()
    target_compile_options(io PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()
