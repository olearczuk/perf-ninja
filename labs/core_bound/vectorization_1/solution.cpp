#include "solution.hpp"

#include <algorithm>
#include <cassert>
#include <type_traits>

namespace simd {
using score_t = int16_t;
using simd_score_t = std::array<score_t, sequence_count_v>;
using simd_sequence_t = std::array<simd_score_t, sequence_size_v>;

// Tranposing results in more cache-friendly layout in memory when processing
// sequences 'concurrently' (via SIMD).
simd_sequence_t transpose(const std::vector<sequence_t> &sequences) {
  assert(sequences.size() == sequence_count_v);
  simd_sequence_t transposed{};

  for (auto i = 0; i < sequences.size(); ++i) {
    for (int j = 0; j < sequences[i].size(); ++j) {
      assert(sequences[i].size() == sequence_size_v);
      transposed[j][i] = sequences[i][j];
    }
  }
  return transposed;
}

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment_simd(std::vector<sequence_t> const &sequences1,
                                std::vector<sequence_t> const &sequences2) {
  result_t result{};

  using column_t = std::array<simd_score_t, sequence_size_v + 1>;

  const auto transposed1 = transpose(sequences1);
  const auto transposed2 = transpose(sequences2);

  /*
   * Initialise score values.
   */
  const score_t gap_open{-11};
  const score_t gap_extension{-1};
  const score_t match{6};
  const score_t mismatch{-4};

  /*
   * Setup the matrix.
   * Note we can compute the entire matrix with just one column in memory,
   * since we are only interested in the last value of the last column in the
   * score matrix.
   */
  column_t score_column{};
  column_t horizontal_gap_column{};
  simd_score_t last_vertical_gap{};

  /*
   * Initialise the first column of the matrix.
   */
  for (int i = 0; i < sequence_count_v; ++i) {
    horizontal_gap_column[0][i] = gap_open;
    last_vertical_gap[i] = gap_open;
  }

  for (size_t i = 1; i < score_column.size(); ++i) {
    for (int k = 0; k < sequence_count_v; ++k) {
      score_column[i][k] = last_vertical_gap[k];
      horizontal_gap_column[i][k] = last_vertical_gap[k] + gap_open;
      last_vertical_gap[k] += gap_extension;
    }
  }

  /*
   * Compute the main recursion to fill the matrix.
   */
  for (unsigned col = 1; col <= transposed2.size(); ++col) {
    simd_score_t last_diagonal_score{};
    for (int k = 0; k < sequence_count_v; ++k) {
      last_diagonal_score[k] = score_column[0][k];  // Cache last diagonal score
                                                    // to compute this cell.
      score_column[0][k] = horizontal_gap_column[0][k];
      last_vertical_gap[k] = horizontal_gap_column[0][k] + gap_open;
      horizontal_gap_column[0][k] += gap_extension;
    }

    for (unsigned row = 1; row <= transposed1.size(); ++row) {
      simd_score_t best_cell_score{};
      for (int k = 0; k < sequence_count_v; ++k) {
        // Compute next score from diagonal direction with
        // match/mismatch.
        best_cell_score[k] =
            last_diagonal_score[k] +
            (transposed1[row - 1][k] == transposed2[col - 1][k] ? match
                                                                : mismatch);
      }

      for (int k = 0; k < sequence_count_v; ++k) {
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score[k] = std::max(best_cell_score[k], last_vertical_gap[k]);
        best_cell_score[k] =
            std::max(best_cell_score[k], horizontal_gap_column[row][k]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_score[k] = score_column[row][k];
        score_column[row][k] = best_cell_score[k];
        // Compute the next values for vertical and horizontal gap.
        best_cell_score[k] += gap_open;
        last_vertical_gap[k] += gap_extension;
        horizontal_gap_column[row][k] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gap[k] =
            std::max(last_vertical_gap[k], best_cell_score[k]);
        horizontal_gap_column[row][k] =
            std::max(horizontal_gap_column[row][k], best_cell_score[k]);
      }
    }
  }

  for (int k = 0; k < sequence_count_v; ++k) {
    result[k] = score_column.back()[k];
  }

  return result;
}
}  // namespace simd

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment_seq(std::vector<sequence_t> const &sequences1,
                               std::vector<sequence_t> const &sequences2) {
  result_t result{};

  for (size_t sequence_idx = 0; sequence_idx < sequences1.size();
       ++sequence_idx) {
    using score_t = int16_t;
    using column_t = std::array<score_t, sequence_size_v + 1>;

    sequence_t const &sequence1 = sequences1[sequence_idx];
    sequence_t const &sequence2 = sequences2[sequence_idx];

    /*
     * Initialise score values.
     */
    score_t gap_open{-11};
    score_t gap_extension{-1};
    score_t match{6};
    score_t mismatch{-4};

    /*
     * Setup the matrix.
     * Note we can compute the entire matrix with just one column in memory,
     * since we are only interested in the last value of the last column in the
     * score matrix.
     */
    column_t score_column{};
    column_t horizontal_gap_column{};
    score_t last_vertical_gap{};

    /*
     * Initialise the first column of the matrix.
     */
    horizontal_gap_column[0] = gap_open;
    last_vertical_gap = gap_open;

    for (size_t i = 1; i < score_column.size(); ++i) {
      score_column[i] = last_vertical_gap;
      horizontal_gap_column[i] = last_vertical_gap + gap_open;
      last_vertical_gap += gap_extension;
    }

    /*
     * Compute the main recursion to fill the matrix.
     */
    for (unsigned col = 1; col <= sequence2.size(); ++col) {
      score_t last_diagonal_score =
          score_column[0];  // Cache last diagonal score to compute this cell.
      score_column[0] = horizontal_gap_column[0];
      last_vertical_gap = horizontal_gap_column[0] + gap_open;
      horizontal_gap_column[0] += gap_extension;

      for (unsigned row = 1; row <= sequence1.size(); ++row) {
        // Compute next score from diagonal direction with match/mismatch.
        score_t best_cell_score =
            last_diagonal_score +
            (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score = std::max(best_cell_score, last_vertical_gap);
        best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_score = score_column[row];
        score_column[row] = best_cell_score;
        // Compute the next values for vertical and horizontal gap.
        best_cell_score += gap_open;
        last_vertical_gap += gap_extension;
        horizontal_gap_column[row] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
        horizontal_gap_column[row] =
            std::max(horizontal_gap_column[row], best_cell_score);
      }
    }

    // Report the best score.
    result[sequence_idx] = score_column.back();
  }

  return result;
}

result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
#ifdef SOLUTION
  return simd::compute_alignment_simd(sequences1, sequences2);
#else
  return compute_alignment_seq(sequences1, sequences2);
#endif
}