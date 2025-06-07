This is a lab about [loop interchange](https://en.wikipedia.org/wiki/Loop_interchange).

[Matrix multiplication](https://en.wikipedia.org/wiki/Matrix_multiplication) is an important building block for many numerical algorithms. In this lab assignment, we compute the integer power of a given real square matrix.
The binary representation of the power significantly reduces the number of matrix operations. Still, the code has a major performance flaw. Your job is to find it out.

## Solution
`multiply` function is the hot path of computing the power.
```c++
void multiply(Matrix &result, const Matrix &a, const Matrix &b) {
  zero(result);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        result[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}
```
You can see that `i` is always present at the first dimension so it should be the outermost-loop for optimal memory access patterns.</br>
However, in `b[k][j]`, `k` is used as the first dimension while `k` is in the innermost-loop.<br/>
Switching `j` and `k` loop is going to result in more efficient access pattern (continuous accesses for each `b[k]`).

## Benchmark
Loop interchange results in ~90% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8929         -0.8930     586092109      62776014     585531442      62645163
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8938         -0.8939     584571866      62080877     583981532      61961341
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8936         -0.8936     591778414      62971334     590883462      62848366
```