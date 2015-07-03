cmake_minimum_required (VERSION 2.6)
project(UTestUtils)
set(CMAKE_CXX_COMPILER c++)
add_executable(utest utest.cpp keys_values.cpp timer.cpp)

include_directories($ENV{HOME}/local/gtest-1.7.0/include/)
target_link_libraries(utest $ENV{HOME}/local/gtest-1.7.0/make/gtest_main.a)
