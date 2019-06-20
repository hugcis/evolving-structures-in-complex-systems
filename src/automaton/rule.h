#include <inttypes.h>

#ifndef RULE_H
#define RULE_H

void build_rule_from_args(uint64_t grule_size,
                          uint8_t rule_array[grule_size],
                          char rule_buf[grule_size + 1],
                          char*, int);

void symmetrize_rule(uint64_t grule_size,
                     uint8_t rule_array[grule_size], int, int);

#endif /* RULE_H */
