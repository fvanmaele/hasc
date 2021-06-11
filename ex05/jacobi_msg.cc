#include <vector>
#include <cassert>
#include <span>
#include <cstdlib>
#include <iomanip>
#include <random>

#include "jacobi_seq.hh"
#include "MessageSystem.hh"
#include "time_experiment.hh"

constexpr int seed = 42;

void thread(const int rank, const int P, int n, int iterations, int inexact,
            std::vector<double>& residuals, std::vector<double>& elapsed)
{
  register_thread(rank);

  int n_local = (n - 2) / P;
  assert((n - 2) == n_local * P);

  // local chunks, includes padding (radius = 1)
  int Ny = n_local + 2;
  int Nx = n;
  int N_local = Ny * Nx;

  std::vector<double> uold_local(N_local, 0.0);
  std::vector<double> unew_local(N_local, 0.0);

  // assign local values (consistent over multiple processes)
  std::mt19937_64 rgen(seed);
  std::uniform_real_distribution<double> dist1(0.0, 1.0);
  rgen.discard(rank * (n - 2) * n_local);

  int offset = 0;
  for (int y = 0; y < Ny; ++y) {
    for (int x = 0; x < Nx; ++x, ++offset) {
      if (y >= 1 && y < Ny - 1 && x >= 1 && x < Nx - 1) {
        uold_local[offset] = unew_local[offset] = dist1(rgen);
      }
    }
  }

  // perform iteration
  auto start = get_time_stamp();
  for (int step = 1; step <= iterations; ++step) {
    // communicate values
    if (P > 1) {
      if (inexact > 0 && step % inexact != 0) {
        goto step;
      }
      ms->barrier(); // ensure data is available on each thread

      // first interior row (send)
      std::vector<double> lsend(Nx);
      for (int i = 0, k = Nx; k < 2*Nx; ++i, ++k) {
        lsend[i] = uold_local[k];
      }

      // last interior row (send)
      std::vector<double> rsend(Nx);
      for (int i = 0, k = N_local - 2*Nx; k < N_local - Nx; ++i, ++k) {
        rsend[i] = uold_local[k];
      }

      // lower boundary row (recv)
      std::vector<double> lrecv(Nx, 0.0);

      // upper boundary column (recv)
      std::vector<double> rrecv(Nx, 0.0);

      // two-sided communication (send/recv order)
      int rnbr = rank+1;
      int lnbr = rank-1;

      if (rank == 0) {
        // right neighbor exchange (S/R)
        send(1, rsend);
        recv(1, rrecv);
      }
      else if (rank == P - 1) {
        // left neighbor exchange (R/S)
        recv(P - 2, lrecv);
        send(P - 2, lsend);
      }
      else if (rank % 2 == 0) {
        // left neighbor exchange (S/R)
        send(lnbr, lsend);
        recv(lnbr, lrecv);

        // right neighbor exchange (S/R)
        send(rnbr, rsend);
        recv(rnbr, rrecv);
      }
      else {
        // right neighbor exchange (R/S)
        recv(rnbr, rrecv);
        send(rnbr, rsend);

        // left neighbor exchange (R/S)
        recv(lnbr, lrecv);
        send(lnbr, lsend);
      }

      // assemble recieved values (no array recv/send available)
      for (int i = 0; i < Nx; ++i) {
        uold_local[i] = lrecv[i]; // left boundary column
      }
      for (int i = N_local - Nx, k = 0; i < N_local; ++i, ++k) {
        uold_local[i] = rrecv[k]; // right boundary column
      }
    } else {
      // nothing to do
    }

  step:
    // print inner regions
#ifdef PRINT_STEPS
    for (int k = 0; k < P; k++) {
      if (rank == k) {
        std::cout << "RANK: " << rank << std::endl;
        std::cout << "STEP: " << step << std::endl;
        std::cout << "GHOST (L): " << std::endl;
        for (int i = 0; i < Nx; ++i) {
          std::cout << std::fixed << uold_local[i] << " ";
        }
        std::cout << "\nINTERIOR: " << std::endl;
        for (int i = Nx, k = 0; i < N_local - Nx; ++i, ++k) {
          if (k > 0 && k % Nx == 0)
            std::cout << std::endl;
          std::cout << std::fixed << uold_local[i] << " ";
        }
        std::cout << "\nGHOST (R): " << std::endl;
        for (int i = N_local - Nx; i < N_local; ++i) {
          std::cout << std::fixed << uold_local[i] << " ";
        }
        std::cout << "\n----\n";
      }
      ms->barrier();
    }
#endif

    // compute current step
    for (int y = 1; y < Ny - 1; y++) {
      for (int x = 1; x < Nx - 1; x++) {
        int index = y*Nx + x;
        unew_local[index] = 0.25 * (uold_local[index - Nx] + uold_local[index - 1] +
                                    uold_local[index + 1] + uold_local[index + Nx]);
      }
    }
    using std::swap;
    swap(uold_local, unew_local);
  }

  // measure vector updates
  auto stop = get_time_stamp();
  elapsed.at(rank) = get_duration_seconds(start, stop);

#ifdef PRINT_RESULT
  for (int k = 0; k < P; k++) {
    if (rank == k) {
      for (int i = Nx, k = 0; i < N_local - Nx; ++i, ++k) {
        if (k > 0 && k % Nx == 0)
          std::cout << std::endl;
        std::cout << std::fixed << uold_local[i] << " ";
      }
      std::cout << std::endl;
    }
    ms->barrier();
  }
#endif

  // compute local residual
  residuals.at(rank) = jacobi_residual(n, uold_local.data(), 1, Ny-1, 1, Nx-1);
}

int main(int argc, char** argv)
{
  int P = 4;
  int n = 1024;
  int iterations = 1000;
  int inexact = 0;

  if (argc < 5) {
    std::cerr << "jacobi_msg <threads> <n> <iterations> <inexact>" << std::endl;
    std::exit(1);
  }
  P = std::strtol(argv[1], nullptr, 10);
  n = std::strtol(argv[2], nullptr, 10);
  iterations = std::strtol(argv[3], nullptr, 10);
  inexact = std::strtol(argv[4], nullptr, 10);

  if (P < 1) {
    std::cerr << "at least one thread required" << std::endl;
    std::exit(1);
  }
  initialize_message_system(P); // assign ms
  std::cerr << "P = " << P << std::endl;

  std::vector<double> residuals(P);
  std::vector<double> elapsed_time(P);
  std::vector<std::thread> threads;
  std::cout << "N,";
  std::cout << "inexact_steps,";
  std::cout << "message_passing_" << P;
  std::cout << std::endl;

  for (int i = 0; i < P; ++i) {
    threads.push_back(std::thread{thread,i, P, n, iterations, inexact,
                                  std::ref(residuals), std::ref(elapsed_time)});
  }
  for (int i = 0; i < P; ++i) {
    threads[i].join();
  }
  std::cout << n * n;
  std::cout << "," << inexact;
  double updates = iterations;
  updates *= (n - 2) * (n - 2);
  double elapsed = elapsed_time[0];
  std::cout << "," << updates / elapsed / 1e9;
  std::cout << std::endl;

  std::cerr << "residual (local): " << std::endl;
  for (int i = 0; i < P; ++i) {
    std::cerr << i << ": " << residuals[i] << std::endl;
  }
  std::cerr << "residual: " << std::endl;
  double res = 0;
  for (int i = 0; i < P; ++i) {
    res += residuals[i]*residuals[i];
  }
  res = std::sqrt(res);
  std::cerr << res << std::endl;
}
