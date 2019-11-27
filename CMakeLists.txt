project(TSISA_LAB_06)
cmake_minimum_required (VERSION 3.8)

add_executable (TSISA_LAB_06 "main.cpp")

target_include_directories(
TSISA_LAB_06
PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)
