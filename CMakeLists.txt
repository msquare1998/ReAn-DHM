# ============================================================================================
#  Re-An method for measuring lnZ of the 1D dimerzied Heisenberg model
#  Reference: https://arxiv.org/abs/2403.08642
#  @Yiming_Ding, Westlake University
#  Last updated: Mar 19, 2024
# ============================================================================================
cmake_minimum_required(VERSION 3.17)
project(thread_DHM)
set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3)
set(THE_SOURCE "./src/main.cpp")
set(OUTPUT_EXE "run")
include_directories("./include")
add_executable(${OUTPUT_EXE} ${THE_SOURCE})
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
target_link_libraries(${OUTPUT_EXE} ${CMAKE_THREAD_LIBS_INIT})