#include <random>

#include "solution.h"

void init(std::vector<S> &arr) {
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(minRandom, maxRandom - 1);

  for (int i = 0; i < N; i++) {
    int random_int1 = distribution(generator);
    int random_int2 = distribution(generator);

    arr[i] = S(random_int1, random_int2);
  }
}
