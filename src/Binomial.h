#pragma once

constexpr uint64_t binomial(int n, int k) {
  if (k > n - k) {
    k = n - k;
  }
  uint64_t result = 1;
  for (int i = 0; i < k; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }
  return result;
}
