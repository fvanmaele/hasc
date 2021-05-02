#include <iostream>
#include <vector>
#include <string>
#include "time_experiment.hh"
#include "transpose.hh"

// initialize square matrix
void initialize(int n, double *A)
{
  for (int i = 0; i < n * n; i++)
    A[i] = i;
}

// conduct experiment X with transposeX
class Experiment1
{
  int n;
  double *A, *B;

public:
  // construct an experiment
  Experiment1(int n_) : n(n_)
  {
    std::cout << "Exp1: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp1: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp1: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment1()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose1(n, A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// conduct experiment X with transposeX
class Experiment2
{
  int n;
  double *A, *B;

public:
  // construct an experiment
  Experiment2(int n_) : n(n_)
  {
    std::cout << "Exp2: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp2: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp2: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment2()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose2(n, A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// conduct experiment X with transposeX
template <int M, int N>
class Experiment3
{
  int n;
  double *A, *B;

public:
  // construct an experiment
  Experiment3(int n_) : n(n_)
  {
    std::cout << "Exp3: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp3: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp3: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment3()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose3<M, N>(n, A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// conduct experiment X with transposeX
template <int M, int N>
class Experiment4
{
  int n;
  double *A, *B;

public:
  // construct an experiment
  Experiment4(int n_) : n(n_)
  {
    std::cout << "Exp4: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp4: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp4: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment4()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose4<M, N>(n, A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// conduct experiment X with transposeX
template <int M>
class Experiment5
{
  int n, b;
  double *A, *B;

public:
  // construct an experiment
  Experiment5(int n_, int b_) : n(n_), b(b_)
  {
    std::cout << "Exp5: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp5: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp5: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment5()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose5<M>(n, std::min(n, b), A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// conduct experiment X with transposeX
template <int M>
class Experiment6
{
  int n, b;
  double *A, *B;

public:
  // construct an experiment
  Experiment6(int n_, int b_) : n(n_), b(b_)
  {
    std::cout << "Exp6: " << n << std::endl;
    A = new double[n * n];
    B = new double[n * n];
    initialize(n, A);
    initialize(n, B);
    if (((size_t)A) % 64 != 0)
    {
      std::cout << "Exp6: A not aligned to 64 " << std::endl;
    }
    if (((size_t)B) % 64 != 0)
    {
      std::cout << "Exp6: B not aligned to 64 " << std::endl;
    }
  }
  ~Experiment6()
  {
    delete[] A;
    delete[] B;
  }
  // run an experiment; can be called several times
  void run() const
  {
    transpose6<M>(n, std::min(n, b), A, B);
  }
  // report number of operations for one run
  double operations() const
  {
    return n * n;
  }
};

// main function runs the experiments and outputs results as csv
int main(int argc, char **argv)
{
  std::vector<int> sizes; // vector with problem sizes to try
  // for (int i=16; i<=16384; i*=2) sizes.push_back(i);
  for (int i = 24; i <= 25000; i *= 2)
    sizes.push_back(i);

  std::vector<std::string> expnames; // name of experiment

  // experiment 1
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
  // experiment 2
  expnames.push_back("vanilla-strided-write");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth2;
  for (auto n : sizes)
  {
    Experiment2 e(n);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth2.push_back(result);
    std::cout << result << std::endl;
  }
  // // experiment 3
  // expnames.push_back("block-consecutive-write-24x4");
  // std::cout << expnames.back() << std::endl;
  // std::vector<double> bandwidth3;
  // for (auto n : sizes)
  //   {
  //     Experiment3<24,4> e(n);
  //     auto d = time_experiment(e,1000000);
  //     double result = d.first*e.operations()*2*sizeof(double)/d.second*1e6/1e9;
  //     bandwidth3.push_back(result);
  //     std::cout << result << std::endl;
  //   }
  // // experiment 4
  // expnames.push_back("block-strided-write-4x24");
  // std::cout << expnames.back() << std::endl;
  // std::vector<double> bandwidth4;
  // for (auto n : sizes)
  //   {
  //     Experiment4<24,4> e(n);
  //     auto d = time_experiment(e,1000000);
  //     double result = d.first*e.operations()*2*sizeof(double)/d.second*1e6/1e9;
  //     bandwidth4.push_back(result);
  //     std::cout << result << std::endl;
  //   }
  // experiment 5
  expnames.push_back("block-consecutive_write-24x4");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5;
  for (auto n : sizes)
  {
    Experiment5<4> e(n, 24);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-48x4");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5b;
  for (auto n : sizes)
  {
    Experiment5<4> e(n, 48);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5b.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-96x4");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5c;
  for (auto n : sizes)
  {
    Experiment5<4> e(n, 48);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5c.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-rowx4");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5d;
  for (auto n : sizes)
  {
    Experiment5<4> e(n, n);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5d.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-24x8");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5e;
  for (auto n : sizes)
  {
    Experiment5<8> e(n, 24);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5e.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-48x8");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5f;
  for (auto n : sizes)
  {
    Experiment5<8> e(n, 48);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5f.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-96x8");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5g;
  for (auto n : sizes)
  {
    Experiment5<8> e(n, 96);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5g.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-rowx8");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5h;
  for (auto n : sizes)
  {
    Experiment5<8> e(n, n);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5h.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-24x24");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5i;
  for (auto n : sizes)
  {
    Experiment5<24> e(n, 24);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5i.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-48x24");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5j;
  for (auto n : sizes)
  {
    Experiment5<24> e(n, 48);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5j.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-96x24");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5k;
  for (auto n : sizes)
  {
    Experiment5<24> e(n, 96);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5k.push_back(result);
    std::cout << result << std::endl;
  }

  expnames.push_back("block-consecutive_write-rowx24");
  std::cout << expnames.back() << std::endl;
  std::vector<double> bandwidth5l;
  for (auto n : sizes)
  {
    Experiment5<24> e(n, n);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth5l.push_back(result);
    std::cout << result << std::endl;
  }

  // experiment 6
  std::cout << expnames.back() << std::endl;
  expnames.push_back("block-strided-writed-4x24");
  std::vector<double> bandwidth6;
  for (auto n : sizes)
  {
    Experiment6<4> e(n, 24);
    auto d = time_experiment(e, 1000000);
    double result = d.first * e.operations() * 2 * sizeof(double) / d.second * 1e6 / 1e9;
    bandwidth6.push_back(result);
    std::cout << result << std::endl;
  }

  // output results
  // Note: size of TLB mentioned in https://www.realworldtech.com/haswell-cpu/5/
  std::cout << "N";
  for (std::string s : expnames)
    std::cout << ", " << s;
  std::cout << std::endl;
  for (int i = 0; i < sizes.size(); i++)
  {
    std::cout << sizes[i];
    std::cout << ", " << bandwidth1[i];
    std::cout << ", " << bandwidth2[i];
    // std::cout << ", " << bandwidth3[i];
    // std::cout << ", " << bandwidth4[i];
    std::cout << ", " << bandwidth5[i];
    std::cout << ", " << bandwidth5b[i];
    std::cout << ", " << bandwidth5c[i];
    std::cout << ", " << bandwidth5d[i];
    std::cout << ", " << bandwidth5e[i];
    std::cout << ", " << bandwidth5f[i];
    std::cout << ", " << bandwidth5g[i];
    std::cout << ", " << bandwidth5h[i];
    std::cout << ", " << bandwidth5i[i];
    std::cout << ", " << bandwidth5j[i];
    std::cout << ", " << bandwidth5k[i];
    std::cout << ", " << bandwidth5l[i];
    std::cout << ", " << bandwidth6[i];
    std::cout << std::endl;
  }

  return 0;
}
