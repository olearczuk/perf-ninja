#include "solution.hpp"

#include <algorithm>
#include <array>
#include <iostream>

unsigned getSumOfDigits(unsigned n) {
  unsigned sum = 0;
  while (n != 0) {
    sum = sum + n % 10;
    n = n / 10;
  }
  return sum;
}

// Task: lookup all the values from l2 in l1.
// For every found value, find the sum of its digits.
// Return the sum of all digits in every found number.
// Both lists have no duplicates and elements placed in *random* order.
// Do NOT sort any of the lists. Do NOT store elements in a hash_map/sets.

// Hint: Traversing a linked list is a long data dependency chain:
//       to get the node N+1 you need to retrieve the node N first.
//       Think how you can execute multiple dependency chains in parallel.
unsigned solution_old(List *l1, List *l2) {
  unsigned retVal = 0;

  List *head2 = l2;
  // O(N^2) algorithm:
  while (l1) {
    unsigned v = l1->value;
    l2 = head2;
    while (l2) {
      if (l2->value == v) {
        retVal += getSumOfDigits(v);
        break;
      }
      l2 = l2->next;
    }
    l1 = l1->next;
  }

  return retVal;
}

template <int ArraySize>
unsigned solution_new(List *l1, List *l2) {
  unsigned retVal = 0;

  std::array<unsigned, ArraySize> arr{};

  while (l1 != nullptr) {
    // store next ArraySize values from l1
    int array_size = 0;
    while (l1 != nullptr && array_size < ArraySize) {
      arr[array_size++] = l1->value;
      l1 = l1->next;
    }

    // iterate over l2 and check whether elements exist in the array
    const List *l = l2;
    while (l != nullptr) {
      auto end_it = arr.cbegin() + array_size;
      auto val = l->value;
      if (auto it = std::find(arr.cbegin(), end_it, val); it != end_it) {
        retVal += getSumOfDigits(val);
      }
      l = l->next;
    }
  }

  return retVal;
}

unsigned solution(List *l1, List *l2) {
#ifdef SOLUTION
  return solution_new<8>(l1, l2);
#else
  return solution_old(l1, l2);
#endif
}
