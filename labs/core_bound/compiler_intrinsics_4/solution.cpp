
#include "solution.h"

#include <immintrin.h>

#include <cassert>
#include <cmath>

#include "const.h"

std::vector<short> mandelbrotOld(int image_width, int image_height) {
  const auto data_width = image_width + 2;
  const auto data_height = image_height + 2;
  const auto diameter_y = kDiameterX / image_width * image_height;
  const auto min_x = kCenterX - kDiameterX / 2;
  const auto max_x = kCenterX + kDiameterX / 2;
  const auto min_y = kCenterY - diameter_y / 2;
  const auto max_y = kCenterY + diameter_y / 2;
  std::vector<short> result(data_width * data_height);
  auto result_idx = 0;
  for (auto py = 0; py < data_height; ++py) {
    for (auto px = 0; px < data_width; ++px) {
      const auto c_x = std::lerp(min_x, max_x, 1.0 * px / data_width);
      const auto c_y = std::lerp(min_y, max_y, 1.0 * py / data_height);
      auto z_x = 0.0;
      auto z_y = 0.0;
      auto iter_cnt = 0;
      for (; iter_cnt < kMaxIterations; ++iter_cnt) {
        const auto z_xx = z_x * z_x;
        const auto z_yy = z_y * z_y;
        if (z_xx + z_yy > kSquareBound) {
          break;
        }
        const auto z_xy = z_x * z_y;
        z_x = z_xx - z_yy + c_x;
        z_y = z_xy + z_xy + c_y;
      }
      result[result_idx++] = iter_cnt;
    }
  }
  return result;
}

std::vector<short> mandelbrotAVX2(int image_width, int image_height) {
  const auto data_width = image_width + 2;
  const auto data_height = image_height + 2;
  const auto diameter_y = kDiameterX / image_width * image_height;
  const auto min_x = kCenterX - kDiameterX / 2;
  const auto max_x = kCenterX + kDiameterX / 2;
  const auto min_y = kCenterY - diameter_y / 2;
  const auto max_y = kCenterY + diameter_y / 2;
  std::vector<short> result(data_width * data_height, kMaxIterations);
  auto result_idx = 0;

  __m256d kSquareBound_vec = _mm256_set1_pd(kSquareBound);
  for (auto py = 0; py < data_height; ++py) {
    const auto c_y = std::lerp(min_y, max_y, 1.0 * py / data_height);
    __m256d c_y_vec = _mm256_set1_pd(c_y);

    // This loop works in the same way as the old approach, but processes
    // 4 pixels in parallel using AVX2 intrinsics.
    int px = 0;
    for (; px + 4 < data_width; px += 4, result_idx += 4) {
      double c_x[4];
      for (int i = 0; i < 4; ++i) {
        c_x[i] = std::lerp(min_x, max_x, 1.0 * (px + i) / data_width);
      }
      __m256d c_x_vec = _mm256_load_pd(c_x);

      __m256d z_x_vec = _mm256_setzero_pd();
      __m256d z_y_vec = _mm256_setzero_pd();

      unsigned int greater_than_bound_mask = 0;
      for (short i = 0; greater_than_bound_mask != 0x0F && i < kMaxIterations;
           ++i) {
        __m256d z_xx_vec = _mm256_mul_pd(z_x_vec, z_x_vec);
        __m256d z_yy_vec = _mm256_mul_pd(z_y_vec, z_y_vec);

        // Check which pixels are greater than the square bound
        __m256d gt_mask = _mm256_cmp_pd(_mm256_add_pd(z_xx_vec, z_yy_vec),
                                        kSquareBound_vec, _CMP_GT_OQ);

        unsigned int mask = _mm256_movemask_pd(gt_mask);
        // This is necessary to avoid
        // processing pixels that were already greater than the bound in
        // previous iterations.
        unsigned int mask_without_previously_processes =
            mask & ~greater_than_bound_mask;
        greater_than_bound_mask |= mask;
        while (mask_without_previously_processes != 0) {
          // Find the first index where the escape condition is met
          int idx = __builtin_ctz(mask_without_previously_processes);
          result[result_idx + idx] = i;

          // Clear the bit for that index
          mask_without_previously_processes &=
              mask_without_previously_processes - 1;
        }

        // Update z_x and z_y
        __m256d z_xy_vec = _mm256_mul_pd(z_x_vec, z_y_vec);
        z_x_vec = _mm256_add_pd(_mm256_sub_pd(z_xx_vec, z_yy_vec), c_x_vec);
        z_y_vec = _mm256_add_pd(_mm256_add_pd(z_xy_vec, z_xy_vec), c_y_vec);
      }
    }

    // Process remaining pixels in the same way as in old approach
    for (; px < data_width; ++px) {
      const auto c_x = std::lerp(min_x, max_x, 1.0 * px / data_width);
      const auto c_y = std::lerp(min_y, max_y, 1.0 * py / data_height);
      auto z_x = 0.0;
      auto z_y = 0.0;
      auto iter_cnt = 0;
      for (; iter_cnt < kMaxIterations; ++iter_cnt) {
        const auto z_xx = z_x * z_x;
        const auto z_yy = z_y * z_y;
        if (z_xx + z_yy > kSquareBound) {
          break;
        }
        const auto z_xy = z_x * z_y;
        z_x = z_xx - z_yy + c_x;
        z_y = z_xy + z_xy + c_y;
      }
      result[result_idx++] = iter_cnt;
    }
  }
  return result;
}

std::vector<short> mandelbrot(int image_width, int image_height) {
#ifdef SOLUTION
  return mandelbrotAVX2(image_width, image_height);
#else
  return mandelbrotOld(image_width, image_height);
#endif
}