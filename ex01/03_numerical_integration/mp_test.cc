#include <limits>
#include <stdexcept>
#include <string>
#include <cstdio>

#include "functions.hh"
#include "midpoint_rule.hh"

#define REQUIRE(expr) Require(expr, __LINE__)

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
    } else {
        std::string str = "test success at line number " + std::to_string(line);
        std::printf("%s\n", str.c_str());
    }  
}

int main()
{
    double mp0 = midpoint_rule_seq(0, 1, 10, [](double x) { return x; });
    REQUIRE(Approx(mp0, 0.5));

    Prim1 F;
    double ig1 = F(1) - F(0); // comparison value
    Func1 f;
    double mp1 = midpoint_rule_seq(0, 1, 10, f);
    REQUIRE(std::abs(mp1 - ig1) < 0.01);
    mp1 = midpoint_rule_seq(0, 1, 100, f);
    REQUIRE(std::abs(mp1 - ig1) < 0.001);

    Prim2 G;
    double ig2 = G(1) - G(0);
    Func2 g;
    double mp2 = midpoint_rule_seq(0, 1, 10, g);
    REQUIRE(std::abs(mp2 - ig2) < 0.1);
    mp2 = midpoint_rule_seq(0, 1, 100, g);
    REQUIRE(std::abs(mp2 - ig2) < 0.01);

    // XXX: add some tests for vectorized implementation
}