cmake_minimum_required(VERSION 3.5.0)
project(CMathUtils VERSION 0.1.0 LANGUAGES C)

# 指定包含目录
include_directories(include)

# 查找 src 目录下的所有 .c 文件
file(GLOB SRC_FILES src/*.c)
# file(GLOB HEAD_FILES include/*.h)

# 添加库
add_library(CMathUtils ${SRC_FILES})  # 添加所有找到的源文件

# 添加可执行文件
add_executable(test_cmu_matrix demo/test_cmu_matrix.c)
target_link_libraries(test_cmu_matrix PRIVATE CMathUtils)

# 添加可执行文件
add_executable(test_cmu_rotation demo/test_cmu_rotation.c)
target_link_libraries(test_cmu_rotation PRIVATE CMathUtils)