cmake_minimum_required(VERSION 3.16)
project(mcpricer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)

add_library(mcpricer INTERFACE)
target_include_directories(mcpricer INTERFACE include)
target_link_libraries(mcpricer INTERFACE OpenMP::OpenMP_CXX)
if(MSVC)
  target_compile_options(mcpricer INTERFACE /W4 /permissive-)
  target_compile_definitions(mcpricer INTERFACE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(mcpricer INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_executable(price src/price.cpp)
target_link_libraries(price PRIVATE mcpricer)

add_executable(bench bench/bench.cpp)
target_link_libraries(bench PRIVATE mcpricer)

add_executable(tests tests/tests.cpp)
target_link_libraries(tests PRIVATE mcpricer)

enable_testing()
add_test(NAME tests COMMAND tests)
