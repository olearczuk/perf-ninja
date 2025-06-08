Loop tiling (blocking) is an important technique that you can use to speed up code that is working with multi-dimensional arrays. If one of the memory access patterns on your array is column-wise, or if in the code you are accessing the same data several times in the loop, this technique can be very beneficial for the performance. It is often seen in matrix multiplication and matrix rotation operations, to speed them up.

Every time the CPU loads a new element of a matrix, it also fetches a few neighboring elements (cache line) belonging to the same row. If matrices are big and you are accessing a matrix column-wise, performance of your code may suffer from poor cache utilization. Because by the time you access the second element in the first row, it's no longer in the cache since it was replaced by the cache lines with elements from other rows of the matrix.

So, instead of going through the whole matrix at once, you can split it into small chunks, which entirely fit into a CPU cache. By processing matrix in blocks (tiles), you are reusing the elements of the matrix which are in the CPU cache and this will give your code a speed boost. Picking the right value for the TILE_SIZE is experimental and depends both on the HW architecture and the algorithm itself. Hint: you can use Roofline Performance analysis (in Intel Advisor or other tools) to determine what's limiting performance of the loop.

Authored-by: @ibogosavljevic

## Solution
The problem is a relatively simple program.<br/>
```c++
using MatrixOfDoubles = std::vector<std::vector<double>>;

bool solution(MatrixOfDoubles &in, MatrixOfDoubles &out) {
  int size = in.size();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      out[i][j] = in[j][i];
    }
  }
  return out[0][size - 1];
}
```
Accesses to `out[i][j]` and `in[j][i]` exhibit poor spatial locality. Each inner loop iteration jumps between rows, which leads to inefficient use of cache lines—especially for `in[j][i]`, where j varies rapidly. This causes frequent cache misses, as adjacent memory loaded with each cache line is rarely reused.
However, this loop does not take advantage of this by jumping around in memory (e.g. fetching entirely different `in[j]` in each inner loop iteration).</br>
This can be solved by [loop tiling technique](https://en.wikipedia.org/wiki/Loop_nest_optimization).<br/>
```c++
bool solution(MatrixOfDoubles &in, MatrixOfDoubles &out) {
  static constexpr int CACHE_LINE_SIZE = 64;
  static constexpr int TILE_SIZE = 8;

  static_assert(sizeof(double) == 8);
  static_assert(sizeof(std::vector<double>) == 24);
  // Ensure that TILE_SIZE is a multiple of CACHE_LINE_SIZE in both dimensions.
  // This way we ensure each tile fits int exactly X cache lines.
  static_assert((TILE_SIZE * sizeof(double)) % CACHE_LINE_SIZE == 0);
  // TILE_SIZE * 24 = 192 bytes per tile row in this dimension.
  // This does not ensure cache alignment, but gives an idea of total size
  // touched per tile block.
  static_assert((TILE_SIZE * sizeof(std::vector<double>)) % CACHE_LINE_SIZE ==
                0);

  int size = in.size();
  for (int i = 0; i < size; i += TILE_SIZE) {
    for (int j = 0; j < size; j += TILE_SIZE) {
      for (int ii = i; ii < std::min(size, i + TILE_SIZE); ++ii) {
        for (int jj = j; jj < std::min(size, j + TILE_SIZE); ++jj) {
          out[ii][jj] = in[jj][ii];
        }
      }
    }
  }
  return out[0][size - 1];
}
```
As you can see, now we're processing both arrays in chunks which take advantage of cache locality (for both vectors `in` and `out`).


## Benchmark
Loop tiling optimization results in ~50% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.4994         -0.5001          7999          4004          7985          3991
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.4953         -0.4952          7944          4009          7918          3997
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.5026         -0.5032          7891          3925          7874          3912
```