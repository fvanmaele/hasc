#include "time_experiment.hh"
#include <vector>
#include <iostream>
#include <random>

static double sum = 0.0;

// Set up experiments to measure the following things:

class Exp1 {
    int n;
public:
    Exp1(int n_) 
        : n(n_) {}

    // 1. An empty loop that iterates from 0 to n-1 without doing anything.
    void run() const {
        for (int i = 0; i < n; ++i);
    }
};

class Exp2 {
    // precision of vector is lower than reduction value for better numerical stability
    const std::vector<float>& v;
public:
    Exp2(const std::vector<float>& v_)
        : v(v_) {}

    // 2. A loop that accumulates something in a global variable.
    void run() const {
        sum = 0; // reset reduction value
        for (auto&& c : v)
            sum += c;
    }
};

class Exp3 {
    // type of the vector is unspecified; use the type of the filled values
    std::vector<size_t>& v;
public:
    Exp3(std::vector<size_t>& v_)
        : v(v_) {}

    // 3. A loop that fills a vector of size n with the integers 0, ..., n-1;
    void run() const {
        for (size_t i = 0; i < v.size(); ++i)
            v[i] = i;
    }
};

int main()
{
    int N = 2;
    int N_lim = 2 << 19;
    int seed = 42;
    std::mt19937_64 gen(seed);
    std::normal_distribution<float> d{ 10, 2 };

    printf("Results for Experiment 1:\n");
    for (int i = N; i <= N_lim; i *= 2) {
        Exp1 f(i);
        
        auto te = time_experiment<Exp1>(f);
        double time_per_run = te.second != 0 ? te.first / te.second : 0;
        std::cout << "N: " << i << ", repetitions: " << te.second << ", runtime(ms): " << time_per_run << "\n";
    }
    
    printf("\nResults for Experiment 2:\n");
    for (int i = N; i <= N_lim; i *= 2) {
        std::vector<float> v(i);
        for (size_t i = 0; i < v.size(); ++i)
            v[i] = d(gen);
        Exp2 f(v);
        
        auto te = time_experiment<Exp2>(f);
        double time_per_run = te.second != 0 ? te.first / te.second : 0;
        std::cout << "N: " << v.size() << ", repetitions: " << te.second << ", runtime(ms): " << time_per_run << "\n";
    }

    printf("\nResults for Experiment 3:\n");
    for (int i = N; i <= N_lim; i *= 2) {
        // note: filled with 0 on initialization
        std::vector<size_t> v(i);
        Exp3 f(v);

        auto te = time_experiment<Exp3>(f);
        double time_per_run = te.second != 0 ? te.first / te.second : 0;
        std::cout << "N: " << v.size() << ", repetitions: " << te.second << ", runtime(ms): " << time_per_run << "\n";
    }

    return 0;
}