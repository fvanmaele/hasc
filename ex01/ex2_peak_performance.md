# Exercise 2 - Peak Performance

## (a) Peak performance

- According to Intel specifications, the Xeon E3-1270 V2 has a maximum bandwidth of **25.6 GBytes/s**. In this section we have run the STREAM benchmark and try to approach this theoretical value.
- According to [APP Metrics for Intel Microprocessors][xeon-app] (Adjusted Peak Performance), the Xeon E3-1270 V2 has an operation throughput of **112 GFLOPS**.

In the sections below, we try to achieve similar numbers using the STREAM and LINPACK benchmarks.

### Memory bandwidth - STREAM

----
Note: we follow *The McCalpin STREAM benchmark: How to do it right and interpret the results* in this section.

----

The SLURM benchmark is compiled with `GCC 10.2` and the following compile options:
```
-Ofast -ffreestanding -fopenmp -fargument-noalias -march=native
```

as well as `-DSTREAM_TYPE=float` and `-DSTream_ARRAY_SIZE=80000000`. This resulted in the following output:
```
Function    Best Rate MB/s  Avg time     Min time     Max time
Copy:           13020.8     0.049986     0.049152     0.051135
Scale:          12998.1     0.049866     0.049238     0.051032
Add:            14552.5     0.067150     0.065968     0.069288
Triad:          14508.3     0.067043     0.066169     0.068346
```
Reducing the amount of threads to 4 and pinning them to physical cores (`OMP_NUM_THREADS=4 OMP_PLACES=4 OMPresulted in only a minor improvement:
```
Function    Best Rate MB/s  Avg time     Min time     Max time
Copy:           13274.7     0.049042     0.048212     0.050534
Scale:          13198.3     0.049362     0.048491     0.050989
Add:            14739.8     0.065763     0.065130     0.066468
Triad:          14810.7     0.065962     0.064818     0.067510
```
A single core configuration (`OMP_NUM_THREADS=1`) only resulted in a small decrease in bandwidth:
```
Function    Best Rate MB/s  Avg time     Min time     Max time
Copy:           13039.1     0.049355     0.049083     0.049872
Scale:          12694.9     0.051131     0.050414     0.051835
Add:            13027.2     0.074467     0.073692     0.075813
Triad:          13074.8     0.074531     0.073424     0.075440
```

If `-ffreestanding` is disabled, *Copy* uses `memcpy()` which in turn uses *streaming store* instructions. This results in a higher repoted bandwidth. See the section below.

### Memory bandwidth - Cache strategy

Because Sandy Bridge/Ivy Bridge use a write-allocate cache strategy, these numbers above need to be multiplied by 1.5 (for Copy and Scale) and 1.33 (for Add and Triad):

|Operation|Adjusted bandwidth|
| ---- | ---- |
| Copy | 19912.05 MB/s |
| Scale | 19797.45 MB/s |
| Add | 19603.9 MB/s |
| Triad | 19698.2 MB/s |

These adjustments are not required with streaming store instructions (direct stores from registers to memory), e.g. as supported by the Intel compiler:
```
Function    Best Rate MB/s  Avg time     Min time     Max time
Copy:           18742.2     0.069278     0.068295     0.071303
Scale:          18684.8     0.069072     0.068505     0.069587
Add:            19500.1     0.099350     0.098461     0.100135
Triad:          19414.3     0.099882     0.098896     0.101120
```

### Operation throughput - LINPACK

The LINPACK benchmark is included with Intel MKL. For 4 cores and 4 threads, we have for solving linear equations
```
=================== Timing linear equation system solver ===================

Size   LDA    Align. Time(s)    GFlops   Residual     Residual(norm) Check
1000   1000   4      0.014      48.6866  1.257328e-12 3.754416e-02   pass
1000   1000   4      0.012      57.6993  1.257328e-12 3.754416e-02   pass
1000   1000   4      0.011      62.2827  1.257328e-12 3.754416e-02   pass
1000   1000   4      0.014      48.6477  1.257328e-12 3.754416e-02   pass
2000   2000   4      0.081      65.9905  4.359901e-12 3.364650e-02   pass
2000   2000   4      0.088      60.6247  4.359901e-12 3.364650e-02   pass
5000   5008   4      1.112      74.9659  2.456835e-11 3.270604e-02   pass
5000   5008   4      1.100      75.7835  2.456835e-11 3.270604e-02   pass
10000  10000  4      8.440      79.0153  9.005856e-11 3.034498e-02   pass
10000  10000  4      9.423      70.7711  9.005856e-11 3.034498e-02   pass
15000  15000  4      32.367     69.5282  2.043429e-10 3.095934e-02   pass
15000  15000  4      30.399     74.0306  2.043429e-10 3.095934e-02   pass
18000  18008  4      54.244     71.6885  2.980030e-10 3.151309e-02   pass
18000  18008  4      48.011     80.9942  2.980030e-10 3.151309e-02   pass
20000  20016  4      65.886     80.9604  4.040749e-10 3.456156e-02   pass
20000  20016  4      66.305     80.4487  4.040749e-10 3.456156e-02   pass
22000  22008  4      89.076     79.7030  4.772430e-10 3.394123e-02   pass
22000  22008  4      89.733     79.1193  4.772430e-10 3.394123e-02   pass
25000  25000  4      141.623    73.5608  5.844116e-10 3.221971e-02   pass
25000  25000  4      150.963    69.0099  5.844116e-10 3.221971e-02   pass
26000  26000  4      168.973    69.3526  6.229394e-10 3.183351e-02   pass
```

Thus ~81 Gflops in this particular benchmark, which is ~72% of the theoretical GFlops.


## (b) Approximation

### Memory bandwidth

The Xeon E3-1270 v2 supports dual-channel DDR3 memory, with a bus width of 64 bits (8 bytes) and bus frequency of 1600 MHz. The memory bandwidth is thus computed as:
```
1600 MHz (frequency) * 2 (dual-channel) * 8 (bus width) = 25600 MB/s
```

### Operation throughput

According to [wikichip][xeon-wikichip]:

> Overall, Sandy Bridge achieves double the FLOP of Nehalem, capable of performing 256-bit multiply + 256-bit add each cycle. 
> That is, Sandy Bridge can sustain 16 single-precision FLOP/cycle or 8 double-precision FLOP/cycle. 

This results in an operation throughput of:
```
3.50 GHz (frequency) * 8 (D-FLOP/cycle) * 4 (cores) = 112 GFLOPS
```

[xeon-app]: https://www.intel.com/content/dam/support/us/en/documents/processors/APP-for-Intel-Xeon-Processors.pdf
[xeon-wikichip]: https://en.wikichip.org/wiki/intel/microarchitectures/sandy_bridge_(client)
[intel-oneapi]: https://software.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top.html