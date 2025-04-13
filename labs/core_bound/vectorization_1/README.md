[Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment) is an important algorithm in many bioinformatics applications and pipelines. The goal of the alignment is to gain insights about their biological  relation. In particular, one is interested how the sequences diverged from a common ancestor by evolutionary events like point mutations or insertions and deletions in the respective sequences.
This problem, however, has quadratic complexity and optimizing it can have a great benefit in many applications.
Since many bioinformatic problems start with the alignment of millions of short sequence pieces of length 150 to 300 symbols, we can gain great performance improvements by using SIMD vectors. In this lab you will learn how the algorithm can be improved by transforming the data layout and exposing SIMD computations.

## Solution
The original `compute_alignment_seq` sequentially iterates over provided `sequences1` and `sequences2` (top-level loop) and performs computations.<br/>
Each sequence is computed independently which means we can take advantage of SIMD instructions. <br/>
This can be achieved by expanding each value within a loop to `N` values where `N` is number of sequences. <br/>
This way there is no data dependency between loop iterations. <br/>
Additionally, transposing the sequences results in more cache-friendly layout.

## Benchmark
This change results in ~50x speedup.
```bash
$ python ../../../tools/check_speedup.py -lab_path $(pwd) -num_runs 3 -v
Benchmark                              Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------
BM_compute_alignment                -0.9793         -0.9793       2733236         56568       2730485         56519
Benchmark                              Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------
BM_compute_alignment                -0.9814         -0.9814       2732671         50930       2729913         50876
Benchmark                              Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------
BM_compute_alignment                -0.9814         -0.9814       2728265         50845       2725507         50803
```