cmake_minimum_required(VERSION 3.2)

project(common)

add_library(common common.cpp)

set_target_properties(common PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(common
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
)

if (MSVC)
    target_compile_options(common PRIVATE
        /W4
    )
else ()
    target_compile_options(common PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()
