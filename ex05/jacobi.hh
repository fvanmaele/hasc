#ifndef GLOBAL_CONTEXT_HH
#define GLOBAL_CONTEXT_HH
#include "jacobi_seq.hh"
#include "jacobi_vcl.hh"
#include "jacobi_par.hh"
#include <memory>

constexpr double TOL = 1e-6;
constexpr int TOL_CHECK = 40;

struct GlobalContext {
  // input data
  int n;          // nxn lattice of points
  int iterations; // number of iterations to do
  double *u0;     // the initial guess
  double *u1;     // temporary vector

  // output data
  GlobalContext(int n_) : n(n_) {}
};

std::pair<int, double> jacobi_vanilla(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_vanilla_kernel(context->n, context->iterations,
                                         uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

std::pair<int, double> jacobi_vectorized(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_vectorized_kernel(context->n, context->iterations,
                                            uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

template <int B>
std::pair<int, double> jacobi_blocked(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_blocked_kernel<B>(context->n, context->iterations,
                                            uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

template <int K>
std::pair<int, double> jacobi_wave(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_wave_kernel<K>(context->n, context->iterations,
                                         uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

template <int K>
std::pair<int, double> jacobi_vectorized_wave(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_vectorized_wave_kernel<K>(context->n, context->iterations,
                                                    uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

std::pair<int, double> jacobi_parallel(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  auto [it, res] = jacobi_parallel_kernel(context->n, context->iterations,
                                          uold, unew, TOL, TOL_CHECK);
  // result should be in u1
  if (context->u1 != uold) {
    std::swap(context->u0, context->u1);
  }
  return std::make_pair(it, res);
}

#endif // GLOBAL_CONTEXT_HH
