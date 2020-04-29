cmake_minimum_required(VERSION 3.2)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(common)
add_subdirectory(graph)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

project(blackholes LANGUAGES CXX)

add_library(blackholes blackholes.cpp)

target_link_libraries(blackholes PUBLIC
    algorithm
    common
    graph
    io
    PRIVATE "${OpenMP_CXX_FLAGS}"
)

target_compile_options(blackholes PRIVATE "${OpenMP_CXX_FLAGS}")

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
