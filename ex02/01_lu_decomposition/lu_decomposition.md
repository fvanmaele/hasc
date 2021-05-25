# Exercise 01 - LU decomposition and Computational Intensity

## Computational intensity

**Q:  What is the computational intensity of the LU decomposition of $A$ if you use standard Gauß elimination without blocking?**

A:  Computing an LU decomposition with standard Gauss elimination requires $2/3 n^3$ floating-point operations, ignoring lower-order terms. (See [Wikipedia][lu-decomposition])

* The machine intensity is for $P_m = 112$ GFlops\sec and $M_m = 25.6$ GBytes\sec equal to:
  
    $I_m = P_m/M_m = 112 / 25.6 = 4.375  \text{  (flops\textbackslash byte)}$

* By definition of the computational intensity, we have with $P = 2/3 n^3$ flops and $M = n^2 \cdot 8$ bytes loaded (using double-precision arithmetic):

    $I = P/M = (2/3 n^3 ) / (8 n^2) = 1/12 n \text{  (flops\textbackslash byte)}$

  Thus for $n >= 53$, we have $I >= 4.416$ and the algorithm is compute-bound. For lower dimensions, it is memory-bound. $n=53$ requires $144 \cdot 53$ bytes or approximately 7.5 KBytes.

* Roofline analysis:

    $P(I) = \min(P_m, IM_m) = \min(112, 1/12 n \cdot 25.6) \approx \min(112, 2.133n) \text{ Gflops\textbackslash sec}$

## LU decomposition with blocking

Now we want to analyze the LU decomposition with blocking, similar to the matrix multiplication example. Make sure you know what is going on and that you could explain the idea with a quick sketch.

## Computational intensity - blocking

**Q: What is the computational intensity for the blocked version?**


## Band matrix

As a last step we want to consider band matrices with band width $k$. Calculate the computational intensity for this case.


[lu-decomposition]: https://en.wikipedia.org/wiki/LU_decomposition
[band-matrix]: https://en.wikipedia.org/wiki/Band_matrix