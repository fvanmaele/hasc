cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define optimization levels
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
                -O3 -march=native -funroll-loops -ffast-math
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

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(benchmark REQUIRED)
find_package(fmt CONFIG REQUIRED)

# add_executable(scalar_product_v0 scalar_product_v0.cc)

# add_executable(scalar_product_v1 scalar_product_v1.cc)
# target_link_libraries(scalar_product_v1 Threads::Threads)

# add_executable(scalar_product_v2 scalar_product_v2.cc)
# target_link_libraries(scalar_product_v2 Threads::Threads)

# add_executable(scalar_product_v3 scalar_product_v3.cc)
# target_link_libraries(scalar_product_v3 Threads::Threads)

# add_executable(scalar_product_v4 scalar_product_v4.cc)
# target_link_libraries(scalar_product_v4 Threads::Threads)

# add_executable(scalar_product_v5 scalar_product_v5.cc)
# target_link_libraries(scalar_product_v5 Threads::Threads)

# add_executable(matmul_seq_v1 matmul_seq_v1.cc)

# add_executable(pointer_chasing pointer_chasing.cc)

# Exercise 00
add_executable(compiler_flags compiler_flags.cc)

add_executable(compiler_flags_bm compiler_flags_bm.cc)
target_link_libraries(compiler_flags_bm 
        benchmark::benchmark benchmark::benchmark_main
)

add_executable(daxpy daxpy.cc)
target_link_libraries(daxpy 
        Threads::Threads 
        fmt::fmt
)

add_executable(daxpy_bm daxpy_bm.cc)
target_link_libraries(daxpy_bm
        Threads::Threads
        benchmark::benchmark benchmark::benchmark_main
)

add_executable(daxpy_bm_2 daxpy_bm_2.cc)
target_link_libraries(daxpy_bm_2
        Threads::Threads
        benchmark::benchmark benchmark::benchmark_main
)