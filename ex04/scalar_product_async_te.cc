#include <iostream>

#include "scalar_product_async.hh"
#include "time_experiment.hh"

#ifndef _WIN32
template <class T>
void doNotOptimizeAway(T&& datum) {
    asm volatile("" : "+r" (datum));
}
#endif

template <typename Callable>
class Experiment
{
private:
    std::span<NumberType> x;
    std::span<NumberType> y;
    Callable sp;
    size_t n;

public:
    Experiment(std::span<NumberType> x_, std::span<NumberType> y_, Callable sp_) :
        x(x_), y(y_), sp(sp_), n(x.size())
    {
        assert(x.size() == y.size());
    }
    void run() const {
#ifdef _WIN32
        throw std::invalid_argument("not implemented for win32");
#else
        doNotOptimizeAway(sp(x, y));
#endif
    }
    double operations() const {
        return 2.0*n;
    }
};

int main() {
    std::vector<size_t> sizes{};
    const size_t n_max = 8<<26;
    const size_t n_min = 8<<3;
    for (size_t n = n_min; n <= n_max; n*=8) {
        sizes.push_back(n);
    }
    std::vector<NumberType> x(n_max);
    std::vector<NumberType> y(n_max);
    sp_init(x, y, -10, 10);
    std::cout << "N,sp_seq,sp_unseq,sp_par,sp_par_unseq,sp_async,sp_async_unseq,sp_packaged_task,sp_openmp,sp_tbb\n";
    for (auto&& i : sizes) {
        auto xsub = std::span{x.begin(), x.begin()+i};
        auto ysub = std::span{y.begin(), y.begin()+i};

        printf("%zu,", i);
        auto e = Experiment(xsub, ysub, sp_async<4>);
        auto d = time_experiment(e);
        double flops = d.first*e.operations()/d.second*1e6/1e9;
        std::cout << "n=" << i << " took " << d.second << " us for " << d.first << " repetitions"
                  << " " << flops << " Gflops/s"
                  << " " << flops*8 << " GByte/s"
                  << std::endl;
    }
    return 0;
}
