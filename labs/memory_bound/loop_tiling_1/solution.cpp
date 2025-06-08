#include "solution.hpp"

#include <algorithm>

bool solutionNew(MatrixOfDoubles &in, MatrixOfDoubles &out) {
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

bool solution(MatrixOfDoubles &in, MatrixOfDoubles &out) {
#ifdef SOLUTION
  return solutionNew(in, out);
#else
  int size = in.size();
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      out[i][j] = in[j][i];
    }
  }
  return out[0][size - 1];
#endif
}
