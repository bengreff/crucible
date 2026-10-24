/* One thread of fixed floating-point and integer work; prints nothing, the caller times it. */
#include <stdint.h>
#include <stdio.h>
int main(void) {
  double x = 1.0, y = 0.5;
  uint64_t h = 1469598103934665603ULL;
  for (long i = 0; i < 400000000L; ++i) {
    x = x * 1.0000000001 + 1e-12;
    y = y * 0.9999999999 + x * 1e-15;
    h = (h ^ (uint64_t)i) * 1099511628211ULL;
  }
  return (x + y + (double)(h & 1)) > 0 ? 0 : 1;
}
