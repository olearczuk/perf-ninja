#include <vector>

// Assume those constants never change
constexpr int N = 1000000;
constexpr int minRandom = 0;
constexpr int maxRandom = 100;

// FIXME: this data structure can be reduced in size
#ifdef SOLUTION
struct S {
  double d;
  short l;
  short i;
  short s;
  bool b;

  bool operator<(const S &s) const { return this->i < s.i; }

  S(int first_value, int second_value)
      : d(static_cast<double>(first_value) / maxRandom),
        l(static_cast<short>(first_value * second_value)),
        i(static_cast<short>(first_value)),
        s(static_cast<short>(second_value)),
        b(first_value < second_value) {}

  // Default constructor required by std::vector
  S() = default;
};

static_assert(sizeof(S) == 16);
#else
struct S {
  int i;
  long long l;
  short s;
  double d;
  bool b;

  bool operator<(const S &s) const { return this->i < s.i; }

  S(int first_value, int second_value)
      : i(first_value),
        l(static_cast<long long>(first_value) * second_value),
        s(static_cast<short>(second_value)),
        d(static_cast<double>(first_value) / maxRandom),
        b(first_value < second_value) {}

  // Default constructor required by std::vector
  S() = default;
};

static_assert(sizeof(S) == 40);
#endif

void init(std::vector<S> &arr);
void solution(std::vector<S> &arr);
