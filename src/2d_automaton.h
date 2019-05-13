#include <math.h>
#include <inttypes.h>

#define N 256
#define STATES 3
#define HORIZON (1)
#define SIDE ( (2 * HORIZON) + 1 )
#define NEIGH_SIZE ( (SIDE * SIDE) - 1 )

#ifndef TWOD_AUTOMATON_H /* Include guard */
#define TWOD_AUTOMATON_H

unsigned long hash(char*);

void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char*);

void symmetrize_rule(uint64_t grule_size, uint8_t rule_array[grule_size]);

void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[], int, int, int, long);

void update_step_general(uint64_t grule_size, uint8_t base[][N], size_t,
                         uint8_t rule[grule_size], uint8_t**);

void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1]);

#endif // TWOD_AUTOMATON_H
