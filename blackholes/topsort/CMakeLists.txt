cmake_minimUm_required(VERSION 3.2)

project(topsort_cli)

add_executable(topsort_cli main.cpp)

add_subdirectory(lib)
set_target_properties(topsort_cli PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_link_libraries(topsort_cli PRIVATE
    blackholes
)

target_include_directories(topsort_cli
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/io"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/graph"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/common"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib/algorithm"
)

if (MSVC)
    target_compile_options(topsort_cli PRIVATE
        /W4
    )
else ()
    target_compile_options(topsort_cli PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()

install(TARGETS topsort_cli RUNTIME DESTINATION bin)
