#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <span>
#include <cassert>
#include <atomic>

#include "nbody_io.hh"
#include "nbody_generate.hh"
#include "time_experiment.hh"

// basic data type for position, velocity, acceleration
typedef double double3[4]; // pad up for later use with SIMD
// basic data type for collapsed loop indices
typedef double double2[2];

/*const double gamma = 6.674E-11;*/
const double G = 1.0;
const double epsilon2 = 1E-10;

// helper variables for threads
const int P = 4; // XXX: fixed number of threads
int count = 0; // count number of threads that arrived at the barrier

std::vector<int> flag(P, 0); // flag indicating waiting thread
std::mutex mx; // mutex for use with the cvs
std::vector<std::condition_variable> cv(P); // for waiting
std::mutex px; // mutex for protecting cout
std::atomic_flag af = ATOMIC_FLAG_INIT; // atomic flag for spin lock

using index = std::pair<int, int>;
using index_v = std::vector<index>;

void barrier (int i)
{
    std::unique_lock<std::mutex> ul{mx};
    count += 1; // one more
    if (count<P)
    {
        // wait on my cv until all have arrived
        flag[i] = 1; // indicate I am waiting
        cv[i].wait(ul, [i]{ return flag[i]==0; }); // wait
    }
    else
    {
        // I am the last one, lets wake them up
        count = 0; // reset counter for next turn

        for (int j=0; j<P; j++)
            if (flag[j]==1)
            {
                flag[j] = 0; // the event
                cv[j].notify_one(); // wake up
            }
    }
}

auto upper_triag_index(int n)
{
    index_v seq;
    seq.reserve((n*(n-1)/2));

    for (int i=0; i<n; ++i) {
        for (int j=i+1; j<n; ++j) {
            seq.push_back({i, j});
        }
    }
    return seq;
}

/** \brief compute acceleration vector from position and masses
 *
 * Executes \sum_{i=0}^{n-1} (n-i-1)*26 = n(n-1)*13
 * flops including 1 division and one square root
 */
void acceleration (std::span<const index> idx, const double3 x[], int n, const double m[], double3 a[])
{
    // local view on a (avoid locks inside a loop)
    std::vector<std::array<double, 4> > a_local(n);
    //double3* a_local = static_cast<double3*>(calloc(n, sizeof(double3))); // zero-initialized

    for (size_t k = 0; k < idx.size(); ++k)
    {
        int i = idx[k].first;
        int j = idx[k].second;
        double d0 = x[j][0]-x[i][0];
        double d1 = x[j][1]-x[i][1];
        double d2 = x[j][2]-x[i][2];

        double r2 = d0*d0 + d1*d1 + d2*d2 + epsilon2;
        double r = sqrt(r2);
        double invfact = G/(r*r2);
        double factori = m[i]*invfact;
        double factorj = m[j]*invfact;

        // Store results as temporaries to reduce contention (critial section)
        //std::unique_lock<std::mutex> lx{px};
        a_local[i][0] += factorj*d0;
        a_local[i][1] += factorj*d1;
        a_local[i][2] += factorj*d2;
        a_local[j][0] -= factori*d0;
        a_local[j][1] -= factori*d1;
        a_local[j][2] -= factori*d2;
    }
    // Update all indices that were updated in this thread
    // Indices that were not updated are zero-additions (calloc)
    std::unique_lock<std::mutex> lx{px};
    for (int i = 0; i < n; ++i)
    {
        a[i][0] += a_local[i][0];
        a[i][1] += a_local[i][1];
        a[i][2] += a_local[i][2];
    }
//    free(a_local);
}

/** \brief do one time step with leapfrog
 *
 * does n*(n-1)*13 + 12n flops
 */
void leapfrog (int rank, const index_v& idx, int n, double dt,
               double3 x[], double3 v[], const double m[], double3 a[])
{
    const int i_begin = n*rank/P;
    const int i_end = n*(rank+1)/P;
    const int N = idx.size();
    const int N_begin = N*rank/P;
    const int N_end = N*(rank+1)/P;

    // Check subdivision bounds
    assert(i_begin * P == n*rank);
    assert(i_end * P = n*(rank+1));
    assert(N_begin * P == N*rank);
    assert(N_end * P == N*(rank+1));

    // update position: 6n flops
    // subdivided evenly between threads
    for (int i = i_begin; i < i_end; i++)
    {
        x[i][0] += dt*v[i][0];
        x[i][1] += dt*v[i][1];
        x[i][2] += dt*v[i][2];
    }

    // save and clear acceleration
    // let each thread assign its partition
    for (int i = i_begin; i < i_end; i++)
    {
        a[i][0] = a[i][1] = a[i][2] = 0.0;
    }
    // required because the subdivision in acceleration is not divided evenly
    // between i, but between i and j (unlike the above)
    barrier(rank);

    // compute new acceleration: n*(n-1)*13 flops
    // subdivision of index space {(i, j) : j>i}
    std::span<const index> idx_local(idx.begin()+N_begin, idx.begin()+N_end);
    acceleration(idx_local, x, n, m, a);

    // update velocity: 6n flops
    // subdivided evenly between threads
    barrier(rank);
    for (int i = i_begin; i < i_end; i++)
    {
        v[i][0] += dt*a[i][0];
        v[i][1] += dt*a[i][1];
        v[i][2] += dt*a[i][2];
    }
}

// Put all work in a single function to avoid overhead from spawning threads
void do_work(int rank, int timesteps, int n, int k, int mod,
             double t, double dt, double m[],
             double3 x[], double3 v[], double3 a[],
             char* name, char* basename, FILE *file)
{
    // initialize timestep and write first file
    std::cout << "step=" << k << " finalstep=" << timesteps << " time=" << t << " dt=" << dt << std::endl;
    auto start = get_time_stamp();
    auto idx = upper_triag_index(n);

    k += 1;
    for (; k<=timesteps; k++)
    {
        leapfrog(rank, idx, n, dt, x, v, m, a);
        t += dt;
        barrier(rank); // finish writing x, v, a for next timestep

        if (rank == 0 && k%mod==0)
        {
            auto stop = get_time_stamp();
            double elapsed = get_duration_seconds(start,stop);
            double flop = mod*(13.0*n*(n-1.0)+12.0*n);
            printf("%g seconds for %g ops = %g GFLOPS\n",elapsed,flop,flop/elapsed/1E9);

            printf("writing %s_%06d.vtk \n",basename,k/mod);
            sprintf(name,"%s_%06d.vtk",basename,k/mod);
            file = fopen(name,"w");
            write_vtk_file_double(file,n,x,v,m,t,dt);
            fclose(file);

            start = get_time_stamp();
        }
    }
}

int main (int argc, char** argv)
{
    int n;              // number of bodies in the system
    double *m;          // array for maasses
    double3 *x;         // array for positions
    double3 *v;         // array for velocites
    double3 *a;         // array for accelerations
    int timesteps;      // final time step number
    int k;              // time step number
    int mod;            // files are written when k is a multiple of mod
    char basename[256]; // common part of file name
    char name[256];     // filename with number
    FILE *file;         // C style file hande
    double t;           // current time
    double dt;          // time step

    // command line for restarting
    if (argc==5)
    {
        sscanf(argv[1],"%s",&basename);
        sscanf(argv[2],"%d",&k);
        sscanf(argv[3],"%d",&timesteps);
        sscanf(argv[4],"%d",&mod);
    }
    else if (argc==6) // command line for starting with initial condition
    {
        sscanf(argv[1],"%s",&basename);
        sscanf(argv[2],"%d",&n);
        sscanf(argv[3],"%d",&timesteps);
        sscanf(argv[4],"%lg",&dt);
        sscanf(argv[5],"%d",&mod);
    }
    else // invalid command line, print usage
    {
        std::cout << "usage: " << std::endl;
        std::cout << "nbody_vanilla <basename> <load step> <final step> <every>" << std::endl;
        std::cout << "nbody_vanilla <basename> <nbodies> <timesteps> <timestep> <every>" << std::endl;
        return 1;
    }

    // set up computation from file
    if (argc==5)
    {
        sprintf(name,"%s_%06d.vtk",basename,k);
        file = fopen(name,"r");
        if (file==NULL)
        {
            std::cout << "could not open file " << std::string(basename) << " aborting" << std::endl;
            return 1;
        }
        n = get_vtk_numbodies(file);
        rewind(file);
        x = static_cast<double3*>(calloc(n,sizeof(double3)));
        v = static_cast<double3*>(calloc(n,sizeof(double3)));
        m = static_cast<double*>(calloc(n,sizeof(double)));
        read_vtk_file_double(file,n,x,v,m,&t,&dt);
        fclose(file);
        k *= mod; // adjust step number
        std::cout << "loaded " << n << "bodies from file " << std::string(basename) << std::endl;
    }

    // set up computation from initial condition
    if (argc==6)
    {
        x = static_cast<double3*>(calloc(n,sizeof(double3)));
        v = static_cast<double3*>(calloc(n,sizeof(double3)));
        m = static_cast<double*>(calloc(n,sizeof(double)));
        //plummer(n,17,x,v,m);
        two_plummer(n,17,x,v,m);
        //cube(n,17,1.0,100.0,0.1,x,v,m);
        std::cout << "initialized " << n << " bodies" << std::endl;
        k = 0;
        t = 0.0;
        printf("writing %s_%06d.vtk \n",basename,k);
        sprintf(name,"%s_%06d.vtk",basename,k);
        file = fopen(name,"w");
        write_vtk_file_double(file,n,x,v,m,t,dt);
        fclose(file);
    }

    // allocate acceleration vector
    a = static_cast<double3*>(calloc(n,sizeof(double3)));

    // leapfrog (int n, double dt, double3 x[], double3 v[], const double m[], double3 a[])
    // do time steps
    std::vector<std::thread> th;
    th.reserve(P);
    for (int i = 0; i < P; ++i) {
        th.push_back(std::thread{do_work, i, timesteps, n, k, mod, t,
                                 dt, m, x, v, a, name, basename, file});
    }
    for (int i = 0; i < P; ++i) {
        th[i].join();
    }
    return 0;
}
