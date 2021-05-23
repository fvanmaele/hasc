#include <iostream>
#include <string>
#include <vectorclass.h>

#include "nbody_io.hh"
#include "nbody_generate.hh"
#include "time_experiment.hh"

// basic data type for position, velocity, acceleration
typedef double double3[4]; // pad up for later use with SIMD

/*const double gamma = 6.674E-11;*/
const double G = 1.0;
const double epsilon2 = 1E-10;

/** \brief compute acceleration vector from position and masses
 *
 * Executes \sum_{i=0}^{n-1} (n-i-1)*26 = n(n-1)*13
 * flops including 1 division and one square root
 */
void acceleration_avx (int n,  const double3 x[], const double m[], double3 a[])
{
  if (n%4 != 0)
    throw std::logic_error("n must be divisible by 4");

  // TODO: check loop unrolling
  for (int i=0; i<n; i++)
  {
    Vec4d xi;
    xi.load(&x[i][0]);

    for (int j=i+1; j<n; j++)
    {
      Vec4d xj;
      xj.load(&x[j][0]); // 0..3
      Vec4d d = xj - xi;
      Vec4d d_sq = d*d;
      double r2 = horizontal_add(d_sq) + epsilon2;
      double r = sqrt(r2);
      double invfact = G/(r*r2);
      Vec4d factori(m[i]*invfact);
      Vec4d factorj(m[j]*invfact);

      Vec4d ai;
      Vec4d aj;
      ai.load(&a[i][0]);
      aj.load(&a[j][0]);
      aj = mul_add(factori, d, ai);
      ai = mul_add(factorj, d, ai);
      ai.store(&a[i][0]);
      aj.store(&a[j][0]);
    }
  }
}

/** \brief do one time step with leapfrog
 *
 * does n*(n-1)*13 + 12n flops
 */
void leapfrog_avx (int n, double dt, double3 x[], double3 v[],
                   const double m[], double3 a[])
{
  if (n%4 != 0)
    throw std::logic_error("n must be divisible by 4");
  // update position: 6n flops
  for (int i=0; i<n; i+=4)
  {
    Vec4d xi_0, xi_1, xi_2, xi_3;
    xi_0.load(&x[i][0]);
    xi_1.load(&x[i+1][0]);
    xi_2.load(&x[i+2][0]);
    xi_3.load(&x[i+3][0]);

    Vec4d vi_0, vi_1, vi_2, vi_3;
    vi_0.load(&v[i][0]);
    vi_1.load(&v[i+1][0]);
    vi_2.load(&v[i+2][0]);
    vi_3.load(&v[i+3][0]);

    Vec4d dt4(dt);
    Vec4d fma_0 = mul_add(dt4, vi_0, xi_0);
    Vec4d fma_1 = mul_add(dt4, vi_1, xi_1);
    Vec4d fma_2 = mul_add(dt4, vi_2, xi_2);
    Vec4d fma_3 = mul_add(dt4, vi_3, xi_3);
    fma_0.store(&x[i][0]);
    fma_1.store(&x[i+1][0]);
    fma_2.store(&x[i+2][0]);
    fma_3.store(&x[i+3][0]);
  }

  // save and clear acceleration
  for (int i=0; i<n; i++)
  {
    a[i][0] = a[i][1] = a[i][2] = 0.0;
  }

  // compute new acceleration: n*(n-1)*13 flops
  acceleration_avx(n,x,m,a);

  // update velocity: 6n flops
  for (int i=0; i<n; i+=4)
  {
    Vec4d vi_0, vi_1, vi_2, vi_3;
    vi_0.load(&v[i][0]);
    vi_1.load(&v[i+1][0]);
    vi_2.load(&v[i+2][0]);
    vi_3.load(&v[i+3][0]);

    Vec4d ai_0, ai_1, ai_2, ai_3;
    ai_0.load(&a[i][0]);
    ai_1.load(&a[i+1][0]);
    ai_2.load(&a[i+2][0]);
    ai_3.load(&a[i+3][0]);

    Vec4d dt4(dt);
    Vec4d fma_0 = mul_add(dt4, ai_0, vi_0);
    Vec4d fma_1 = mul_add(dt4, ai_1, vi_1);
    Vec4d fma_2 = mul_add(dt4, ai_2, vi_2);
    Vec4d fma_3 = mul_add(dt4, ai_3, vi_3);
    fma_0.store(&v[i][0]);
    fma_1.store(&v[i+1][0]);
    fma_2.store(&v[i+2][0]);
    fma_3.store(&v[i+3][0]);
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

  // initialize timestep and write first file
  std::cout << "step=" << k << " finalstep=" << timesteps << " time=" << t << " dt=" << dt << std::endl;
  auto start = get_time_stamp();

  // do time steps
  k += 1;
  for (; k<=timesteps; k++)
  {
    leapfrog_avx(n,dt,x,v,m,a);
    t += dt;
    if (k%mod==0)
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

  return 0;
}
