cmake_minimum_required (VERSION 3.12)
project(hasc)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
#find_package(benchmark REQUIRED)
#find_package(fmt CONFIG REQUIRED)
find_package(vectorclass REQUIRED)

# Define optimization levels
set(OPT_VEC_REPORT ON)
set(OPT_FAST_MATH ON)
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
        -Wall -Wextra -Wpedantic -march=native -funroll-loops
        )
endif()

#add_compile_definitions(NDEBUG)

# Exercise 03

add_executable(nbody_vanilla nbody_vanilla.cc)
target_compile_options(nbody_vanilla PUBLIC -fno-tree-vectorize)

add_executable(nbody_vcl nbody_vcl.cc)
target_link_libraries(nbody_vcl vectorclass::vectorclass)

add_executable(nbody_threaded nbody_threaded.cc)
target_link_libraries(nbody_threaded Threads::Threads)
