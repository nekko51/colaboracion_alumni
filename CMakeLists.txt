cmake_minimum_required(VERSION 3.10)
project(ColaboracionAlumni C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")

find_package(OpenMP REQUIRED)

set(COMMON_SOURCES
    c-code/chain-operations.c
    c-code/correlations.c
    c-code/metropolis.c
    c-code/parsing.c
    c-code/random.c
    c-code/svm.c
    c-code/utils.c
)

add_executable(main_svm c-code/main-svm.c ${COMMON_SOURCES})
target_link_libraries(main_svm PRIVATE OpenMP::OpenMP_C m)
target_include_directories(main_svm PRIVATE c-code)

add_executable(main_sim c-code/main.c ${COMMON_SOURCES})
target_link_libraries(main_sim PRIVATE OpenMP::OpenMP_C m)
target_include_directories(main_sim PRIVATE c-code)