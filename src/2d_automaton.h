#include <math.h>
#include <inttypes.h>

#define N 256

#ifndef TWOD_AUTOMATON_H /* Include guard */
#define TWOD_AUTOMATON_H

unsigned long hash(char*);

void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char*, int);

void symmetrize_rule(uint64_t grule_size,
                     uint8_t rule_array[grule_size], int, int);

void process_rule(uint64_t grule_size, uint8_t rule[grule_size],
                  char rule_buf[], int, int, int, long, int, int);

void update_step_general(uint64_t grule_size, uint8_t base[][N], size_t,
                         uint8_t rule[grule_size], uint8_t**, int, int);

void generate_general_rule(uint64_t grule_size,
                           uint8_t rule_array[grule_size],
                           char rule_buf[grule_size + 1], int, int);

#endif // TWOD_AUTOMATON_H
