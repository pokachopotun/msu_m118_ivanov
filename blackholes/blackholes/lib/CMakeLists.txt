cmake_minimum_required(VERSION 3.2)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(common)
add_subdirectory(graph)

project(blackholes)

add_library(blackholes blackholes.cpp)

target_link_libraries(blackholes PUBLIC
    algorithm
    common
    graph
    io
)

set_target_properties(blackholes PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(blackholes
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/algorithm"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/common"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/graph"
)

if (MSVC)
    target_compile_options(blackholes PRIVATE
        /W4
    )
else ()
    target_compile_options(blackholes PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()
