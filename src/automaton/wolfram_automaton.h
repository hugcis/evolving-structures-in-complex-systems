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
  enum InitMode init; /**< Initialiation parameter */
  size_t timesteps; /**< Number of timesteps for the automaton */
  enum WriteMode write; /**< Wether to write output or not to stdin */
  int grain; /**< Granularity of the compression */
  int radius; /**< Radius of the rule window */
};

void write_to_file(size_t size, size_t rule_size, uint8_t rule[rule_size],
                   int, struct Options1D*, int states);

unsigned long rule_number(int states, size_t rule_size,
                          uint8_t rule[rule_size]);
