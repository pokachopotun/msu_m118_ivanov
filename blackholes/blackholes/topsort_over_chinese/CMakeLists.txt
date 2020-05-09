cmake_minimUm_required(VERSION 3.2)

project(topsort_over_chinese_cli)

add_executable(topsort_over_chinese_cli main.cpp)

set_target_properties(topsort_over_chinese_cli PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_link_libraries(topsort_over_chinese_cli PRIVATE
    blackholes
)

target_include_directories(topsort_over_chinese_cli
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/io"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/graph"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/common"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/algorithm"
)

if (MSVC)
    target_compile_options(topsort_over_chinese_cli PRIVATE
        /W4
    )
else ()
    target_compile_options(topsort_over_chinese_cli PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()

install(TARGETS topsort_over_chinese_cli RUNTIME DESTINATION bin)
