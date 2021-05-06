#ifndef HASC_EX02_TRANSPOSE_AVX_HH
#define HASC_EX02_TRANSPOSE_AVX_HH
#include <vectorclass.h>
#include <cassert>
#include <iosfwd>

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// Matrix accessed in a block col-major fashion
inline void
transpose1_avx(int n, const double* __restrict A, double* __restrict B)
{
  assert(n % 4 == 0);

  for (int i = 0; i < n; i += 4) {
    Vec4d a0, a1, a2, a3;
    Vec4d b0, b1, b2, b3;

    for (int j = 0; j < n; j += 4) {
      // Load row of A (4 entries) into vector register
      a0.load(&A[j*n + i]);
      a1.load(&A[(j + 1)*n + i]);
      a2.load(&A[(j + 2)*n + i]);
      a3.load(&A[(j + 3)*n + i]);

      // Reorder into columns
      b0 = blend4<0, 4, 2, 6>(a0, a1);
      b1 = blend4<1, 5, 3, 7>(a0, a1);
      b2 = blend4<0, 4, 2, 6>(a2, a3);
      b3 = blend4<1, 5, 3, 7>(a2, a3);
      a0 = blend4<0, 1, 4, 5>(b0, b2);
      a1 = blend4<0, 1, 4, 5>(b1, b3);
      a2 = blend4<2, 3, 6, 7>(b0, b2);
      a3 = blend4<2, 3, 6, 7>(b1, b3);

      // Stored transposed row of A into B
      a0.store(&B[i*n + j]);
      a1.store(&B[(i + 1)*n + j]);
      a2.store(&B[(i + 2)*n + j]);
      a3.store(&B[(i + 3)*n + j]);
    }
  }
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// Matrix accessed in a block row-major fashion
inline void
transpose2_avx(int n, const double* __restrict A, double* __restrict B)
{
  for (int j = 0; j < n; j += 4) {
    Vec4d a0, a1, a2, a3;
    Vec4d b0, b1, b2, b3;

    for (int i = 0; i < n; i += 4) {
      // Load row of A (4 entries) into vector register
      a0.load(&A[j*n + i]);
      a1.load(&A[(j + 1)*n + i]);
      a2.load(&A[(j + 2)*n + i]);
      a3.load(&A[(j + 3)*n + i]);

      // Reorder into columns
      b0 = blend4<0, 4, 2, 6>(a0, a1);
      b1 = blend4<1, 5, 3, 7>(a0, a1);
      b2 = blend4<0, 4, 2, 6>(a2, a3);
      b3 = blend4<1, 5, 3, 7>(a2, a3);
      a0 = blend4<0, 1, 4, 5>(b0, b2);
      a1 = blend4<0, 1, 4, 5>(b1, b3);
      a2 = blend4<2, 3, 6, 7>(b0, b2);
      a3 = blend4<2, 3, 6, 7>(b1, b3);

      // Stored transposed row of A into B
      a0.store(&B[i*n + j]);
      a1.store(&B[(i + 1)*n + j]);
      a2.store(&B[(i + 2)*n + j]);
      a3.store(&B[(i + 3)*n + j]);
    }
  }
}

template<int M> // bxb blocks with M blocking in the strided direction
inline void
transpose5_avx(int n, int b, const double* __restrict A, double* __restrict B)
{
  static_assert(M >= 4);
  if (b % M != 0) {
    std::cout << b << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % b != 0) {
    std::cout << n << " is not a multiple of " << b << std::endl;
    exit(0);
  }
  if (b % 4 != 0) {
    std::cout << b << " is not a multiple of " << 4 << std::endl;
    exit(0);
  }

  // first two loops over the blocks
  for (int I = 0; I < n; I += b) {
    for (int J = 0; J < n; J += b) {
      for (int j = J; j < J + b; j += M) {
        for (int i = 0; i < I + b; i += 4) {
          Vec4d a0, a1, a2, a3;
          Vec4d b0, b1, b2, b3;

          for (int jj = 0; jj < M; jj += 4) {
            int jjj = j + jj;

            // Load row of A (4 entries) into vector register
            a0.load(&A[jjj*n + i]);
            a1.load(&A[(jjj + 1)*n + i]);
            a2.load(&A[(jjj + 2)*n + i]);
            a3.load(&A[(jjj + 3)*n + i]);

            // Reorder into columns
            b0 = blend4<0, 4, 2, 6>(a0, a1);
            b1 = blend4<1, 5, 3, 7>(a0, a1);
            b2 = blend4<0, 4, 2, 6>(a2, a3);
            b3 = blend4<1, 5, 3, 7>(a2, a3);
            a0 = blend4<0, 1, 4, 5>(b0, b2);
            a1 = blend4<0, 1, 4, 5>(b1, b3);
            a2 = blend4<2, 3, 6, 7>(b0, b2);
            a3 = blend4<2, 3, 6, 7>(b1, b3);

            // Stored transposed row of A into B
            a0.store(&B[i*n + jjj]);
            a1.store(&B[(i + 1)*n + jjj]);
            a2.store(&B[(i + 2)*n + jjj]);
            a3.store(&B[(i + 3)*n + jjj]);
          }
        }
      }
    }
  }
}

template<int M> // bxb blocks with M blocking in the strided direction
inline void
transpose6_avx(int n, int b, const double* __restrict A, double* __restrict B)
{
  static_assert(M >= 4);
  if (b % M != 0) {
    std::cout << b << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % b != 0) {
    std::cout << n << " is not a multiple of " << b << std::endl;
    exit(0);
  }
  if (b % 4 != 0) {
    std::cout << b << " is not a multiple of " << 4 << std::endl;
    exit(0);
  }

  // first two loops over the blocks
  for (int I = 0; I < n; I += b) {
    for (int J = 0; J < n; J += b) {
      for (int j = J; j < J + b; j += M) {
        for (int jj = 0; jj < M; jj += 4) {
          Vec4d a0, a1, a2, a3;
          Vec4d b0, b1, b2, b3;

          for (int i = 0; i < I + b; i += 4) {
            int jjj = j + jj;

            // Load row of A (4 entries) into vector register
            a0.load(&A[jjj*n + i]);
            a1.load(&A[(jjj + 1)* n + i]);
            a2.load(&A[(jjj + 2)* n + i]);
            a3.load(&A[(jjj + 3)* n + i]);

            // Reorder into columns
            b0 = blend4<0, 4, 2, 6>(a0, a1);
            b1 = blend4<1, 5, 3, 7>(a0, a1);
            b2 = blend4<0, 4, 2, 6>(a2, a3);
            b3 = blend4<1, 5, 3, 7>(a2, a3);
            a0 = blend4<0, 1, 4, 5>(b0, b2);
            a1 = blend4<0, 1, 4, 5>(b1, b3);
            a2 = blend4<2, 3, 6, 7>(b0, b2);
            a3 = blend4<2, 3, 6, 7>(b1, b3);

            // Stored transposed row of A into B
            a0.store(&B[i*n + jjj]);
            a1.store(&B[(i + 1)*n + jjj]);
            a2.store(&B[(i + 2)*n + jjj]);
            a3.store(&B[(i + 3)*n + jjj]);
          }
        }
      }
    }
  }
}

#endif