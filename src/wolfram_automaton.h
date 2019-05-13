#include <stdlib.h>
#include <stdint.h>

enum InitMode {
               ONE,
               RANDOM,
               RAND_SMALL
};

enum WriteMode {
                NO_WRITE,
                WRITE_STEP
};

struct Options1D {
  enum InitMode init;
  size_t timesteps;
  enum WriteMode write;
};

void write_to_file(size_t size, size_t rule_size, uint8_t rule[rule_size],
                   int, int, struct Options1D*);

unsigned long rule_number(size_t rule_size, uint8_t rule[rule_size]);
