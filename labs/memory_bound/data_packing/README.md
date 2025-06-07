# Data packing

This is a lab about data packing.

You can decrease the memory traffic of the application if you pack the data more efficiently.
Some of the ways to do that include:

* Eliminate compiler-added padding.
* Use types that require less memory or less precision e.g. (int -> short, double -> float).
* Use bitfields to pack the data even further.

## Solution
We start with this struct:
```c++
struct S {
  int i;       // 4 bytes + 4 bytes of padding
  long long l; // 8 bytes
  short s;     // 2 byte + 6 bytes of padding
  double d;    // 8 bytes
  bool b;      // 1 byte
               // extra 7 bytes of padding
};
static_assert(sizeof(S) == 40);
```

### Step 1: Field Reordering
By reorganizing fields to minimize padding, we reduce the struct size by 40%:
```c++
struct S {
  double d;     // 8 bytes
  long long l;  // 8 bytes
  int i;        // 4 bytes
  short s;      // 2 bytes
  bool b;       // 1 byte
                // extra 1 byte of padding
};
static_assert(sizeof(S) == 24);
```

### Step 2: Type Optimization
Analyzing the data usage:
```c++
constexpr int minRandom = 0;
constexpr int maxRandom = 100;

void init(std::vector<S> &arr) {
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(minRandom, maxRandom - 1);

  for (int i = 0; i < N; i++) {
    int random_int1 = distribution(generator);
    int random_int2 = distribution(generator);

    arr[i] = S(random_int1, random_int2);
  }
}

// Constructor shows actual value ranges:
S(int first_value, int second_value)
    : i(first_value), // < 100
      l(static_cast<long long>(first_value) * second_value), // < 10,000
      s(static_cast<short>(second_value)),
      d(static_cast<double>(first_value) / maxRandom),
      b(first_value < second_value) {}
```

Since all values are small (product < `10,000`), we can use `short` (max `32,767`) instead of `int` and `long long`.
Final optimized struct:
```c++
struct S {
  double d;   // 8 bytes
  short l;    // 2 bytes
  short i;    // 2 bytes
  short s;    // 2 bytes
  bool b;     // 1 byte
              // extra 1 byte of padding
};
static_assert(sizeof(S) == 16);
```

## Benchmark
The 60% reduction in memory footprint (from 40 to 16 bytes) results in a corresponding ~60% performance improvement.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.6074         -0.6073            32            12            32            12
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.5961         -0.5960            31            13            31            12
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.6054         -0.6053            32            13            32            13
```