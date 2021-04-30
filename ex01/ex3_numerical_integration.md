# Exercise 3

## (a)

See `functions.hh` and `midpoint_rule.hh`.

## (b)

See `functions.hh` and `midpoint_rule.hh`. Vectorization is done with the `vectorclass` library and `Vec4d` (256-bit vector, double, AVX).

## (c)

Benchmark results:
```
----------------------------------------------------------------------
Benchmark                            Time             CPU   Iterations
----------------------------------------------------------------------
BM_midpoint_f1_seq/10             22.3 ns         22.5 ns     29866667
BM_midpoint_f1_seq/64             29.6 ns         29.8 ns     23578947
BM_midpoint_f1_seq/512             218 ns          220 ns      3200000
BM_midpoint_f1_seq/4096           1722 ns         1726 ns       407273
BM_midpoint_f1_seq/32768         13975 ns        14125 ns        49778
BM_midpoint_f1_seq/262144       125197 ns       119978 ns         5600
BM_midpoint_f1_seq/2097152      870484 ns       878514 ns          747
BM_midpoint_f1_seq/10485760    4533546 ns      4492188 ns          160
```
```
BM_midpoint_f2_seq/10             82.3 ns         82.0 ns      8960000
BM_midpoint_f2_seq/64              560 ns          558 ns      1120000
BM_midpoint_f2_seq/512            4526 ns         4450 ns       154483
BM_midpoint_f2_seq/4096          36512 ns        36901 ns        19478
BM_midpoint_f2_seq/32768        287921 ns       284934 ns         2358
BM_midpoint_f2_seq/262144      2333196 ns      2351589 ns          299
BM_midpoint_f2_seq/2097152    19198208 ns     19425676 ns           37
BM_midpoint_f2_seq/10485760   93462329 ns     93750000 ns            7
```
```
BM_midpoint_f1_vec/10             15.0 ns         15.1 ns     49777778
BM_midpoint_f1_vec/64             48.6 ns         47.6 ns     14451613
BM_midpoint_f1_vec/512             330 ns          330 ns      2133333
BM_midpoint_f1_vec/4096           2543 ns         2490 ns       263529
BM_midpoint_f1_vec/32768         20763 ns        20403 ns        34462
BM_midpoint_f1_vec/262144       164778 ns       164958 ns         4073
BM_midpoint_f1_vec/2097152     1323512 ns      1317771 ns          498
BM_midpoint_f1_vec/10485760    6778696 ns      6835938 ns          112
```
```
BM_midpoint_f2_vec/10             47.4 ns         47.1 ns     14933333
BM_midpoint_f2_vec/64              144 ns          143 ns      4480000
BM_midpoint_f2_vec/512            1070 ns         1050 ns       640000
BM_midpoint_f2_vec/4096           8448 ns         8545 ns        89600
BM_midpoint_f2_vec/32768         68818 ns        68011 ns         8960
BM_midpoint_f2_vec/262144       568645 ns       578125 ns         1000
BM_midpoint_f2_vec/2097152     4318716 ns      4261364 ns          154
BM_midpoint_f2_vec/10485760   21630581 ns     21484375 ns           32
```