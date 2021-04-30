#include <limits>
#include <stdexcept>
#include <string>

#include "functions.hh"
#include "midpoint_rule.hh"
#include <iostream>

using namespace hasc_ex01;

bool Approx(double a, double b)
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

int main()
{
    double mp0 = midpoint_rule_seq(0, 1, 10, [](double x) { return x; });
    Require(Approx(mp0, 0.5), __LINE__);

    Prim1 F;
    double ig1 = F(1) - F(0); // comparison value
    Func1 f;
    double mp1 = midpoint_rule_seq(0, 1, 10, f);
    Require(std::abs(mp1 - ig1) < 0.01, __LINE__);
    mp1 = midpoint_rule_seq(0, 1, 100, f);
    Require(std::abs(mp1 - ig1) < 0.001, __LINE__);

    Prim2 G;
    double ig2 = G(1) - G(0);
    Func2 g;
    double mp2 = midpoint_rule_seq(0, 1, 100, g);
    Require(std::abs(mp2 - ig2) < 0.01, __LINE__);
    mp2 = midpoint_rule_seq(0, 1, 1000, g);
    Require(std::abs(mp2 - ig2) < 0.001, __LINE__);
}