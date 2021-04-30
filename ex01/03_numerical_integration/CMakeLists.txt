cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(benchmark REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(fmt CONFIG REQUIRED)
find_package(vectorclass REQUIRED)

# Define optimization levels
set(OPT_VEC_REPORT ON)

# Vectorization reports
if (OPT_VEC_REPORT)
        if (MSVC)
        add_compile_options(
                /Qvec-report:1
        )
        elseif(CMAKE_COMPILER_IS_GNUCXX)
        add_compile_options(
                -fopt-info-vec-all
        )
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -Rpass=loop-vectorize
        )
        endif()
endif()

# Set optimization flags
if (MSVC)
        add_compile_options(
                /W3 /arch:AVX
        )
elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -Wall -Wextra -Wpedantic -march=ivybridge
        )
endif()

#add_compile_definitions(NDEBUG)

# Exercise 01

add_executable(mp_bench mp_bench.cc)
target_link_libraries(mp_bench PRIVATE 
        vectorclass::vectorclass
        OpenMP::OpenMP_CXX
        benchmark::benchmark benchmark::benchmark_main
)

add_executable(mp_test mp_test.cc)
target_link_libraries(mp_test PRIVATE 
        vectorclass::vectorclass
        OpenMP::OpenMP_CXX
        benchmark::benchmark benchmark::benchmark_main
)
