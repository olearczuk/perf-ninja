#include "solution.h"

#include <random>

void generateObjectsNew(InstanceArray& array) {
  std::default_random_engine generator(0);
  std::uniform_int_distribution<std::uint32_t> distribution(0, 2);

  std::size_t a_count = 0, b_count = 0, c_count = 0;
  for (std::size_t i = 0; i < N; i++) {
    int value = distribution(generator);
    if (value == 0) {
      ++a_count;
    } else if (value == 1) {
      ++b_count;
    } else {
      ++c_count;
    }
  }

  for (std::size_t i = 0; i < a_count; i++) {
    array.push_back(std::make_unique<ClassA>());
  }
  for (std::size_t i = 0; i < b_count; i++) {
    array.push_back(std::make_unique<ClassB>());
  }
  for (std::size_t i = 0; i < c_count; i++) {
    array.push_back(std::make_unique<ClassC>());
  }
}

void generateObjects(InstanceArray& array) {
#ifdef SOLUTION
  generateObjectsNew(array);
#else
  std::default_random_engine generator(0);
  std::uniform_int_distribution<std::uint32_t> distribution(0, 2);

  for (std::size_t i = 0; i < N; i++) {
    int value = distribution(generator);
    if (value == 0) {
      array.push_back(std::make_unique<ClassA>());
    } else if (value == 1) {
      array.push_back(std::make_unique<ClassB>());
    } else {
      array.push_back(std::make_unique<ClassC>());
    }
  }
#endif
}

// Invoke the `handle` method on all instances in `output`
void invoke(InstanceArray& array, std::size_t& data) {
  for (const auto& item : array) {
    item->handle(data);
  }
}
