cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
#find_package(benchmark REQUIRED)
#find_package(fmt CONFIG REQUIRED)
find_package(vectorclass REQUIRED)

# Define optimization levels
# XXX: Whatever you do here, release/etc will take precedence. CMake == Cursed Make
set(OPT_LEVEL_LOW OFF)
set(OPT_LEVEL_MEDIUM OFF)
set(OPT_LEVEL_HIGH ON)
set(OPT_VEC_REPORT OFF)

# Vectorization reports
if (OPT_VEC_REPORT)
        if (MSVC)
        add_compile_options(
                /Qvec-report:1
        )
        elseif(CMAKE_COMPILER_IS_GNUCXX)
        add_compile_options(
                -fopt-info-vec
        )
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -Rpass=loop-vectorize
        )
        endif()
endif()

# Set optimization flags
if (OPT_LEVEL_LOW)
        if (MSVC)
        add_compile_options(
                /Od
        )
        elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -O0
        )
        endif()
elseif(OPT_LEVEL_MEDIUM)
        if (MSVC)
        add_compile_options(
                /O1 /arch:AVX
                )
        elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -O2 -march=native
        )
        endif()
elseif(OPT_LEVEL_HIGH)
        if (MSVC)
        add_compile_options(
                /O2 /Ob3 /arch:AVX /fp:fast
        )
        elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -O3 -march=native -ffast-math
        )
        endif()
endif()

# Set warning levels
if (MSVC)
        add_compile_options(
                /W3
        )
elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(
                -Wall -Wextra -Wpedantic
        )
endif()

#add_compile_definitions(NDEBUG)

# Exercise 01
add_executable(matmul_seq_v1 matmul_seq_v1.cc)
add_executable(matmul_seq_v2 matmul_seq_v2.cc)
target_link_libraries(matmul_seq_v2 PRIVATE vectorclass::vectorclass)

add_executable(pointer_chasing pointer_chasing.cc)


