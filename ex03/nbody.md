# N-Body problem

## Total energy

We checked this with Python on VTK-exported `.csv` files (see `total_energy.py`). Total energy remained roughly the same, as long as the timesteps were chosen sufficiently small (`t<0.001`).

This proved useful to check the correctness of the [threaded](#Threaded) version.

## Blocking

Blocking either resulted in no performance gains or extremely bad performance (~0.01 GFlops). As such we did not further explore this option.

## Vectorization

Manual vectorization resulted in no performance improvements over auto-vectorization (with GCC 10.2 and `AVX` instructions). The disassembly (`objdump -d`) was also similar in both cases.

Note: `horizontal_add` for the sum `d0*d0 + d1*d1 + d2*d2 + epsilon` can also be implemented with `permute4`:

```cpp
Vec4d d_sq1 = permute4<2, 0, 1, 3>(d_sq);
Vec4d d_sq2 = permute4<1, 2, 0, 3>(d_sq);
Vec4d eps(epsilon2);
Vec4d r2 = d_sq + d_sq1 + d_sq2 + eps;
Vec4d r = sqrt(r2);
Vec4d invfact = G/(r*r2);
Vec4d mi(m[i]);
Vec4d mj(m[j]);
Vec4d factori = mi * invfact;
Vec4d factorj = mj * invfact;
```
Performance was similar in both cases.

## Threading

We used the following approach, resulting in a performance gain of `~4x` for 4 threads.

* Acceleration: distribute the index space `{(i,j) : j > i}` for `0 <= i < n` evenly between threads.
* Velocity, position: distribute the index space `{i | 0 <= i < n}` evenly between threads.
* Use a barrier (cf. `barrier.cc`) after computing each time step to handle data dependencies.

This resulted in better load balancing, but requires an additional barrier (cf. `barrier.cc`) after setting `a = 0` (Note the index space `i` is distributed differently for acceleration and velocity/position, respectively).

* To avoid overhead from spawning threads, all computations were moved inside a single function `do_work()`. Files were printed on the thread with rank `0`.

### Critical section

The critical section in each thread is updating the `a` array. To avoid locks inside the main computation loop (at the cost of additional storage space), we gave each thread a "view" on `a`:

```cpp
double3* a_local = static_cast<double3*>(calloc(n, sizeof(double3))); // zero-initialized
```

and updated it in a seperate loop:

```cpp
std::unique_lock<std::mutex> lx{px};
for (int i = 0; i < n; ++i)
{
    a[i][0] += a_local[i][0];
    a[i][1] += a_local[i][0];
    a[i][2] += a_local[i][0];
}
```