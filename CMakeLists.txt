# 设置最低 CMake 版本
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(flight_sim VERSION 0.1 LANGUAGES C)
# 设置MSVC编译编码
#add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/source-charset:utf-8>")

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_FLAGS /source-charset:utf-8)

aux_source_directory(src SOURCES)
list(REMOVE_ITEM SOURCES src/test.c) # 排除测试文件

# 添加可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    include
)

# # 链接库（如 winmm.lib）
# if(WIN32)
#     target_link_libraries(flight_sim winmm)
# endif()
