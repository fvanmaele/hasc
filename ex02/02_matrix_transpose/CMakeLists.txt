cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(benchmark REQUIRED)
#find_package(fmt CONFIG REQUIRED)
find_package(vectorclass REQUIRED)

# Define optimization levels
set(OPT_VEC_REPORT ON)
set(OPT_FAST_MATH OFF)

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

if (OPT_FAST_MATH)
        if (MSVC)
        add_compile_options(
                /fp:fast
        )
        elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -ffast-math
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
                -Wall -Wextra -Wpedantic -march=native
        )
endif()

#add_compile_definitions(NDEBUG)

# Exercise 02

add_executable(transpose transpose.cc)

add_executable(transpose_avx transpose_avx.cc transpose_avx.hh)
target_link_libraries(transpose_avx PRIVATE vectorclass::vectorclass)

add_executable(transpose_trimatrix transpose_trimatrix.cc trimatrix.hh)
add_executable(transpose_test transpose_test.cc)
target_link_libraries(transpose_test PRIVATE vectorclass::vectorclass)