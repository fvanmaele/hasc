cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)

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
add_executable(pointer_chasing pointer_chasing.cc)

