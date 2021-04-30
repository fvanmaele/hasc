#ifndef HASC_EX01_MIDPOINT_RULE
#define HASC_EX01_MIDPOINT_RULE
#include <type_traits>
#include <cassert>
#include <iostream>
#include "functions.hh"

namespace hasc_ex01
{

/// @brief Approximates integral of a function R->R over an interval [a,b]
/// @tparam Callable is a function object with a double argument
/// @param a is beginning of the interval
/// @param b is end of the interval
/// @param n the number of steps (rectangles) to use
/// @param f is the function to integrate
template <typename Callable>
double midpoint_rule_seq(double a, double b, int n, Callable f)
{
    assert(n >= 1);
    static_assert(std::is_invocable<Callable, double>::value, "not f(double");

    double sum = 0;
    for (int k = 1; k <= n; ++k) {
        double m = a + (2.*k - 1) / (2.*n) * (b - a);
        sum += f(m);
    }
    return (b - a) * sum / n;
}

template <typename Callable>
double midpoint_rule_vec(double a, double b, int n, Callable f)
{
    assert(n >= 4);
    static_assert(std::is_invocable<Callable, Vec4d>::value, "not f(Vec4d");

    double m0 = a + (2.*0 - 1) / (2.*n) * (b - a);
    double m1 = a + (2.*1 - 1) / (2.*n) * (b - a);
    double m2 = a + (2.*2 - 1) / (2.*n) * (b - a);
    double m3 = a + (2.*3 - 1) / (2.*n) * (b - a);
    Vec4d sum_v = Vec4d(m0, m1, m2, m3);
    int n_last_block = (n / 4) * 4;

    for (int k = 4; k <= n_last_block; k += 4) {
        m0 = a + (2.*k - 1) / (2.*n) * (b - a);
        m1 = a + (2.*(k+1) - 1) / (2.*n) * (b - a);
        m2 = a + (2.*(k+2) - 1) / (2.*n) * (b - a);
        m3 = a + (2.*(k+3) - 1) / (2.*n) * (b - a);
        sum_v += f(Vec4d(m0, m1, m2, m3));
    }
    // Reduce partial sums
    double psums[4] __attribute__((aligned(32))); // 256-bit vector
    sum_v.store(psums);
    double sum = 0;
    for (int i = 0; i < 4; ++i)
        sum += psums[i];
    
    // Add last block
    // XXX: Assumes overload for scalar function!
    for (int k = n_last_block+1; k <= n; ++k) {
        double m = a + (2.*k - 1) / (2.*n) * (b - a);
        sum += f(m);
    }
    return (b - a) * sum / n;
}

} // namespace hasc_ex01

#endif // HASC_EX01_MIDPOINT_RULE