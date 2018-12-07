cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 11)

project(AP-Tree)

include_directories(include)
aux_source_directory(lib SOURCES_LIB)
add_library(${PROJECT_NAME} STATIC ${SOURCES_LIB})

aux_source_directory(test SOURCES_TEST)
add_executable(${PROJECT_NAME}Test ${SOURCES_TEST})
target_link_libraries(${PROJECT_NAME}Test ${PROJECT_NAME})
