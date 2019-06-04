#include "utils.h"

uint32_t ipow(int base, int exp)
{
  uint32_t result = 1;
  for (;;)
    {
      if (exp & 1)
        result *= base;
      exp >>= 1;
      if (!exp)
        break;
      base *= base;
    }

  return result;
}
