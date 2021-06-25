#ifndef SEIDEL_HH
#define SEIDEL_HH
#include "seidel_seq.hh"
#include "seidel_tbb.hh"
#include <memory>

constexpr double TOL = 1e-6;
constexpr int TOL_CHECK = 40;

struct GlobalContext {
  // input data
  int n;          // nxn lattice of points
  int iterations; // number of iterations to do
  double *u0;     // the initial guess

  // output data
  GlobalContext(int n_) : n(n_) {}
};

std::pair<int, double> seidel_vanilla(std::shared_ptr<GlobalContext> context) {
  double* u = context->u0;
  auto [it, res] = seidel_vanilla_kernel(context->n, context->iterations, u, TOL, TOL_CHECK);

  return std::make_pair(it, res);
}

template <int B>
std::pair<int, double> seidel_blocked(std::shared_ptr<GlobalContext> context) {
  double* u = context->u0;
  auto [it, res] = seidel_blocked_kernel<B>(context->n, context->iterations, u, TOL, TOL_CHECK);

  return std::make_pair(it, res);
}

// std::pair<int, double> seidel_tbb(std::shared_ptr<GlobalContext> context) {
//   double* u = context->u0;
//   auto [it, res] = seidel_tbb_kernel(context->n, context->iterations, u, TOL, TOL_CHECK);

//   return std::make_pair(it, res);
// }

#endif // SEIDEL_HH
