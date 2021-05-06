# Exercise 01 - LU decomposition and Computational Intensity

## Computational intensity

**Q:  What is the computational intensity of the LU decomposition of $A$ if you use standard Gauß elimination without blocking?**

A:  Computing an LU decomposition with standard Gauss elimination requires $2/3 n^3$ floating-point operations, ignoring lower-order terms. (See [Wikipedia][lu-decomposition])

* By definition of the computational intensity, we have with $P = 2/3 n^3$ flops and $M = n^2$ bytes:

    $I = (2/3 n^3 ) / n^2 = 2/3 n \text{ flops\textbackslash byte}$

* The machine intensity is for $P_m = 112$ GFlops\sec and $M_m = 25.6$ GBytes\sec equal to:
  
    $I_m = 112 / 25.6 = 4.375  \text{ flops\textbackslash byte}$

* Roofline analysis:

    $P(I) = \min(P_m, IM_m) = P_m \min(1, I/I_m) \approx 1.5238$

The algorithm is thus a compute-bound algorithm.

## LU decomposition with blocking

Now we want to analyze the LU decomposition with blocking, similar to the matrix multiplication example. Make sure you know what is going on and that you could explain the idea with a quick sketch.

## Computational intensity - blocking

What is the computational intensity for the blocked version?

## Band matrix

As a last step we want to consider band matrices with band width $k$. Calculate the computational intensity for this case.


[lu-decomposition]: https://en.wikipedia.org/wiki/LU_decomposition
[band-matrix]: https://en.wikipedia.org/wiki/Band_matrix