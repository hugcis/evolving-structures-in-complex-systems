#include <math.h>
#include <inttypes.h>

#define N 256
#define STATES 2
#define HORIZON (1)
#define SIDE ( (2 * HORIZON) + 1 )
#define NEIGH_SIZE ( (SIDE * SIDE) - 1 )
#define GRULE_SIZE  ( (int)pow(STATES, NEIGH_SIZE + 1) )

#ifndef TWOD_AUTOMATON_H /* Include guard */
#define TWOD_AUTOMATON_H

unsigned long hash(char*);

void build_rule_from_args(uint8_t rule_array[GRULE_SIZE],
                          char rule_buf[GRULE_SIZE + 1],
                          char** argv);

void symmetrize_rule(uint8_t rule_array[GRULE_SIZE]);

void process_rule(uint8_t rule[GRULE_SIZE], char rule_buf[],
                  int, int, int);

void update_step_general(uint8_t base[][N], size_t,
                         uint8_t rule[GRULE_SIZE], uint8_t**);

#endif // TWOD_AUTOMATON_H
