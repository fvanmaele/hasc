# Exercise 01 - LU decomposition and Computational Intensity

## Computational intensity

**Q:  What is the computational intensity of the LU decomposition of $A$ if you use standard Gauß elimination without blocking?**

A:  Computing an LU decomposition with standard Gauss elimination requires $2/3 n^3$ floating-point operations, ignoring lower-order terms. (See [Wikipedia][lu-decomposition])

* By definition of the computational intensity, we have with $P = 2/3 n^3$ flops and $M = n^2 * 8$ bytes loaded (using double-precision arithmetic):

    $I = (2/3 n^3 ) / (8 * n^2) = 1/12 n \text{ flops\textbackslash byte}$

  Thus for $n >= 12$, the algorithm is compute-bound and otherwise memory-bound. $n=12$ requires $144 * 8$ bytes or approximately 1 KByte.

* The machine intensity is for $P_m = 112$ GFlops\sec and $M_m = 25.6$ GBytes\sec equal to:
  
    $I_m = 112 / 25.6 = 4.375  \text{ flops\textbackslash byte}$

* Roofline analysis:

    $P(I) = \min(P_m, IM_m) = \min(112, 1/12 n \cdot 25.6) \approx \min(112, 2.133n) \text{ Gflops\textbackslash sec}$


## LU decomposition with blocking

Now we want to analyze the LU decomposition with blocking, similar to the matrix multiplication example. Make sure you know what is going on and that you could explain the idea with a quick sketch.

## Computational intensity - blocking

**Q: What is the computational intensity for the blocked version?**

A: In [Computational intensity](#computational-intensity), we have computed the intensity (under perfect conditions) as a function of the dimension $n$. For the blocked version, we use the block size $m$ instead. In particular, *per block* intensity is lower ($m \ll n$) and block size should be $>= 12x12$.

## Band matrix

As a last step we want to consider band matrices with band width $k$. Calculate the computational intensity for this case.


[lu-decomposition]: https://en.wikipedia.org/wiki/LU_decomposition
[band-matrix]: https://en.wikipedia.org/wiki/Band_matrix