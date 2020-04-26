cmake_minimUm_required(VERSION 3.2)

project(chinese_cli)

add_executable(chinese_cli main.cpp)

set_target_properties(chinese_cli PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

target_link_libraries(chinese_cli PRIVATE
    blackholes
)

target_include_directories(chinese_cli
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/io"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/graph"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/common"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../lib/algorithm"
)

if (MSVC)
    target_compile_options(chinese_cli PRIVATE
        /W4
    )
else ()
    target_compile_options(chinese_cli PRIVATE
        -Wall -Wextra -pedantic -Werror
    )
endif()

install(TARGETS chinese_cli RUNTIME DESTINATION bin)
