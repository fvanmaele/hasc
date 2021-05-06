#include <limits>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <iostream>
#include <new>

#include "transpose.hh"
#include "transpose_avx.hh"

#define REQUIRE(expr) Require(expr, __LINE__)

bool ApproxEqual(double a, double b)
{
    return (std::abs(a-b) < std::numeric_limits<double>::epsilon());
}

void Require(bool expr, int line)
{
    if (!expr) {
        std::string str = "test failure at line number " + std::to_string(line);
        throw std::invalid_argument(str);
    }
}

// initialize square matrix
void initialize(int n, double *A)
{
  for (int i = 0; i < n * n; i++)
    A[i] = i;
}

int main()
{
    const int n = 256;
    double* A = new(std::align_val_t{64}) double[n * n];
    double* B = new(std::align_val_t{64}) double[n * n];
    double* B2 = new(std::align_val_t{64}) double[n * n];
    initialize(n, A);
    initialize(n, B);
    initialize(n, B2);

    transpose1(n, A, B);
    transpose1_avx(n, A, B2);
    for (int i = 0; i < n; ++i) {
        REQUIRE(ApproxEqual(B[i], B2[i]));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            REQUIRE(ApproxEqual(A[n*i + j], B[n*j + i]));
            REQUIRE(ApproxEqual(A[n*i + j], B2[n*j + i]));
        }
    }
    
    initialize(n, A);
    initialize(n, B);
    initialize(n, B2);
    transpose2(n, A, B);
    transpose2_avx(n, A, B2);
    for (int i = 0; i < n; ++i) {
        REQUIRE(ApproxEqual(B[i], B2[i]));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            REQUIRE(ApproxEqual(A[n*i + j], B[n*j + i]));
            REQUIRE(ApproxEqual(A[n*i + j], B2[n*j + i]));
        }
    }
    
    initialize(n, A);
    initialize(n, B);
    initialize(n, B2);
    transpose5<4>(n, 16, A, B);
    transpose5_avx<4>(n, 16, A, B2);
    for (int i = 0; i < n; ++i) {
        REQUIRE(ApproxEqual(B[i], B2[i]));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            REQUIRE(ApproxEqual(A[n*i + j], B[n*j + i]));
            REQUIRE(ApproxEqual(A[n*i + j], B2[n*j + i]));
        }
    }
    
    initialize(n, A);
    initialize(n, B);
    initialize(n, B2);
    transpose6<4>(n, 16, A, B);
    transpose6_avx<4>(n, 16, A, B2);
    for (int i = 0; i < n; ++i) {
        REQUIRE(ApproxEqual(B[i], B2[i]));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            REQUIRE(ApproxEqual(A[n*i + j], B[n*j + i]));
            REQUIRE(ApproxEqual(A[n*i + j], B2[n*j + i]));
        }
    }

    std::cout << "All tests successful" << std::endl;
}