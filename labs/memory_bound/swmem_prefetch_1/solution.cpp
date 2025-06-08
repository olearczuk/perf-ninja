#include "solution.hpp"

static int getSumOfDigits(int n) {
  int sum = 0;
  while (n != 0) {
    sum = sum + n % 10;
    n = n / 10;
  }
  return sum;
}

int solutionNew(const hash_map_t *hash_map, const std::vector<int> &lookups) {
  constexpr int PREFETCH_DISTANCE = 16;

  int result = 0;
  const auto process_index = [&](int i) {
    int val = lookups[i];
    if (hash_map->find(val)) {
      result += getSumOfDigits(val);
    }
  };

  // Process elements while prefetching the next bucket.
  int i = 0;
  for (; i < lookups.size() - PREFETCH_DISTANCE; ++i) {
    process_index(i);
    __builtin_prefetch(
        hash_map->find_bucket_ptr(lookups[i + PREFETCH_DISTANCE]));
  }

  // Process the remaining elements without prefetching
  for (; i < lookups.size(); ++i) {
    process_index(i);
  }

  return result;
}

int solution(const hash_map_t *hash_map, const std::vector<int> &lookups) {
#ifdef SOLUTION
  return solutionNew(hash_map, lookups);
#else
  int result = 0;
  for (int val : lookups) {
    if (hash_map->find(val)) {
      result += getSumOfDigits(val);
    }
  }
  return result;
#endif
}