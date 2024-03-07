#include "Utilities.hpp"

int factorial(int n) {
  int output = 1;
  for(int i = 1; i <= n; i++) {
    output *= i;
  }
  return output;
}