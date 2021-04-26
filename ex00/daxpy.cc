#include <thread>
#include <vector>
#include <algorithm> // for min
#include <numeric> // for iota
#include <cassert>
#include <utility> // for ref
#include <fmt/format.h>

#include "time_experiment.hh"

// Serial version for comparison purposes. As floating point operations are
// not associative, the result of the parallel version may differ due to grouping.
void daxpy_serial(const std::vector<double> &x, std::vector<double> &y, double a)
{
    assert(x.size() == y.size());
    const size_t N = x.size();

    for (size_t i = 0; i < N; ++i) {
        y[i] += a*x[i];
    }
}

// As the elements of the result vector are written directly to disjoint blocks,
// we do not need locking or condition variables as with a scalar product
// and its global sum.
//
// Note: y is passed by non-const reference, requiring std::ref when used
// with std::thread. An alternative is to take a pointer instead of a reference
// (e.g. throughs std::vector::data)
void daxpy(const std::vector<double> &x, std::vector<double> &y, double a,
           int rank, int nproc)
{
    assert(rank < nproc);
    assert(x.size() == y.size());
    const size_t N = x.size();
    
    // Disjoint partition for input and result vectors
    for (size_t i = (N*rank)/nproc; i < (N*(rank+1))/nproc; ++i) {
        y[i] += a*x[i];
    }
}

class Experiment
{
private:
    int n, nproc;
    double a;
    const std::vector<double>& x;
    std::vector<double>& y;
    
public:
    Experiment(int n_, int nproc_, double a_, const std::vector<double> &x_, std::vector<double> &y_) :
        n(n_), nproc(nproc_), a(a_), x(x_), y(y_)
    {}
    void run() const
    {
        std::vector<std::thread> threads;
        threads.reserve(nproc);

        for (int rank = 0; rank < nproc-1; ++rank) {
            threads.push_back(std::thread{daxpy, x, std::ref(y), a, rank, nproc});
        }
        // Use the main thread for computations and spawn one thread less
        daxpy(x, y, a, nproc-1, nproc);
        
        for (int rank = 0; rank < nproc-1; ++rank) {
            threads[rank].join();
        }
    }
    // report number of (floating-point) operations
    double operations() const
    {
        // a*x[i] + y[i]  ->  N-1 additions, N multiplications
        return 2.0*n - 1;
    }
};

int main()
{
    std::vector<std::thread> threads;
    // We choose the amount of threads depending on problem size to reduce
    // overhead on small vectors.
    const int min_size = 32;
    const int hardware_threads = std::thread::hardware_concurrency();
    const int min_threads = hardware_threads != 0 ? hardware_threads : 2;
    const int Na = 32;
    const int Nb = 2 << 20;
    
    for (int N = Na; N <= Nb; N *= 2) {
        // Fill values of x (1, 2, 3, ..., N)
        std::vector<double> x(N);
        std::iota(x.begin(), x.end(), 1);
        
        // Fill values of y (N, N-1, ..., 1)
        // Note: y will have different values over multiple runs
        std::vector<double> y(N);
        std::iota(y.rbegin(), y.rend(), 1);
        
        // Run experiment
        int nproc = std::min(min_threads, (N+min_size-1) / min_size);
        Experiment e(N, nproc, 1.5, x, y);
        
        auto d = time_experiment(e);
        double flops = d.first*e.operations()/d.second*1e6/1e9;
        double Gflops = flops*sizeof(double);

        fmt::print("N: {}, threads: {}, repetitions: {}, total time: {}us, {} Gflops/s, {} GByte/s\n",
                    N, nproc, d.first, d.second, flops, Gflops);
    }
}
