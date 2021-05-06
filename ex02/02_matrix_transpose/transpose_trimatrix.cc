
#include <vector>
#include <iostream>
#include <string>
#include "trimatrix.hh"
#include "time_experiment.hh"

// conduct experiment X with transposeX
class Experiment1
{
    int n;
    TriMatrix<double> *A, *B;

public:
    // construct an experiment
    Experiment1(int n_) : n(n_)
    {
        A = new TriMatrix<double>(n_);
        B = new TriMatrix<double>(n_);

        int k = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                (*A)(i, j) = k++;
                (*B)(i, j) = 0;
            }
        }

        std::cout << "Exp1: " << n << std::endl;
    }
    ~Experiment1()
    {
        delete A;
        delete B;
    }
    // run an experiment; can be called several times
    void run() const
    {
        double* B_lower = B->lower();
        double* B_upper = B->upper();
        double* A_lower = A->lower();
        double* A_upper = A->upper();

        for (ptrdiff_t i = 0, t = A->t(); i < t; ++i)
        {
            B_lower[i] = A_upper[i];
            B_upper[i] = A_lower[i];
        }
    }
    // report number of operations for one run
    double operations() const
    {
        return n * (n - 1); // set upper and lower triangles
    }
};

int main()
{
    std::vector<int> sizes; // vector with problem sizes to try
    for (int i = 24; i <= 25000; i *= 2)
        sizes.push_back(i);

    std::vector<std::string> expnames; // name of experiment
    expnames.push_back("vanilla-consecutive-write");
    std::cout << expnames.back() << std::endl;
    std::vector<double> bandwidth1;
    for (auto n : sizes)
    {
        Experiment1 e(n);
        auto d = time_experiment(e, 1000000);
        double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
        bandwidth1.push_back(result);
        std::cout << result << std::endl;
    }

    std::cout << "N";
    for (std::string s : expnames)
        std::cout << ", " << s;
    std::cout << std::endl;
    for (int i = 0; i < sizes.size(); i++)
    {
        std::cout << sizes[i];
        std::cout << ", " << bandwidth1[i];
        std::cout << std::endl;
    }

    return 0;
}