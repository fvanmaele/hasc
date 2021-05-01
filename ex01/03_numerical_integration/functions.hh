#ifndef HASC_EX01_FUNCTIONS_HH
#define HASC_EX01_FUNCTIONS_HH
#include <cmath>
#include <vectorclass/vectorclass.h>
#include "pow_int.h"

namespace hasc_ex01
{
// f(x) = x^3 - 2x^2 + 3x - 1
struct Func1
{
    double operator()(double x)
    {
        return x*(x*(x - 2) + 3) - 1;
    }
    Vec4d operator()(Vec4d x)
    {
        return x*(x*(x - 2) + 3) - 1;
    }
};

// f(x) = x^0 + x + ... + x^15
struct Func2
{
    double operator()(double x)
    {
        return (x + 1) * (x*x + 1) * (pow_int(x, 4) + 1) * (pow_int(x, 8) + 1);
    }
    Vec4d operator()(Vec4d x)
    {
        Vec4d ones = Vec4d(1, 1, 1, 1);
        return (x + ones) * (x*x + ones) * (pow_const(x, 4) + 1) * (pow_const(x, 8) + 1);
    }
};

// TODO: unoptimized comparison versions
// F(x) = 1/4 x^4  - 2/3 x^3 + 3/2 x^2 - x
struct Prim1
{
    double operator()(double x)
    {
        return pow_int(x, 4)/4 - 2*pow_int(x, 3)/3 + 3*x*x/2 - x;
    }
    Vec4d operator()(Vec4d x)
    {
        return pow_const(x, 4)/4 - 2*pow_const(x, 3)/3 + 3*x*x/2 - x;
    }
};

// F(x) = x + 1/2 x^2 + 1/3 x^3 + ... + 1/16 x^16
struct Prim2
{
    double operator()(double x)
    {
        double sum = x;
        for (int i = 2; i <= 16; ++i) {
            sum += pow_int(x, i) / i;
        }
        return sum;
    }
    Vec4d operator()(Vec4d x)
    {
        Vec4d ones = Vec4d(1, 1, 1, 1);
        return ones + x + pow_const(x, 2)/2 + pow_const(x, 3)/3 + pow_const(x, 4)/4 +
               pow_const(x, 5)/5 + pow_const(x, 6)/6 + pow_const(x, 7)/7 + pow_const(x, 8)/8 +
               pow_const(x, 9)/9 + pow_const(x, 10)/10 + pow_const(x, 11)/11 + pow_const(x, 12)/12 +
               pow_const(x, 13)/13 + pow_const(x, 14)/14 + pow_const(x, 15)/15;
    }
};

} // namespace hasc_ex01

#endif // HASC_EX01_FUNCTIONS_HH