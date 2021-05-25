cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(benchmark REQUIRED)
find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)
#find_package(fmt CONFIG REQUIRED)
#find_package(vectorclass REQUIRED)

# Define optimization levels
set(OPT_VEC_REPORT OFF)
set(OPT_FAST_MATH OFF)
set(OPT_TSAN OFF)
set(OPT_ASAN OFF)

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

# ThreadSanitizer
if (OPT_TSAN)
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-fsanitize=thread -fno-omit-frame-pointer)
        add_link_options(-fsanitize=thread)
    else()
        message("thread sanitizer not supported for this compiler")
    endif()
endif()

# AddressSanitizer
if (OPT_ASAN)
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address)
    else()
        message("address sanitizer not supported for this compiler")
    endif()
endif()

# Set optimization flags
if (MSVC)
    add_compile_options(
        /W3 /arch:AVX
        )
elseif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(
        -Wall -Wextra -Wpedantic -march=native -fno-omit-frame-pointer
        )
endif()

#add_compile_definitions(NDEBUG)

# Exercise 04

add_library(async INTERFACE scalar_product_async.hh)
target_link_libraries(async INTERFACE Threads::Threads TBB::tbb OpenMP::OpenMP_CXX)

add_executable(tree_barrier tree_barrier.cc Barrier.hh)
target_link_libraries(tree_barrier Threads::Threads)

add_executable(scalar_product_async scalar_product_async.cc)
target_link_libraries(scalar_product_async benchmark::benchmark async)

add_executable(scalar_product_async_te scalar_product_async_te.cc time_experiment.hh)
target_link_libraries(scalar_product_async_te async)

add_executable(bm_minimal_reproducer bm_minimal_reproducer.cc)
target_link_libraries(bm_minimal_reproducer benchmark::benchmark Threads::Threads)