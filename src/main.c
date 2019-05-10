#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "2d_automaton.h"
#include "wolfram_automaton.h"

#define STEPS 1000

const char help[] = "Use with either 2d or 1d as first argument";

int main_2d(int argc, char** argv)
{
  char rule_buf[GRULE_SIZE + 1];
  uint8_t rule_array[GRULE_SIZE];

  /* Write steps for a given rule  */
  if (argc == 3) {
    build_rule_from_args(rule_array, rule_buf, argv[2]);
    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("Rule %s\n", rule_buf);

    process_rule(rule_array, rule_buf, 1, 1, 0, STEPS);
    return 0;
  }


  /* Generate compression plots for many rules */
  time_t t;
  srand((unsigned)time(&t));
  FILE* dic_file;
  char* fname;

  for (int i = 0; i < 1000; ++i) {
    /* generate_totalistic_rule(rule_array, rule_buf); */
    for (int v = 0; v < GRULE_SIZE; v++) {
      /* if (rand() % 10 < 3) { */
        /* rule_array[v] = 1; */
      /* } else { */
        /* rule_array[v] = 0; */
      /* } */
      rule_array[v] = (uint8_t)(rand() % STATES);
    }

    symmetrize_rule(rule_array);
    for (int v = 0; v < GRULE_SIZE; v++) {
      sprintf(&rule_buf[v], "%"PRIu8, rule_array[v]);
    }

    asprintf(&fname, "data_2d_%i/%lu.map", STATES, hash(rule_buf));
    dic_file = fopen(fname, "w+");
    fprintf(dic_file, "%s", rule_buf);
    fclose(dic_file);

    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("%i: Rule %s\n", i, rule_buf);
    fflush(stdout);

    process_rule(rule_array, rule_buf, 0, 1, 0, STEPS);
  }
  printf("\n");
  return 0;
}

int main_1d(int argc, char** argv)
{
  uint32_t rule_size = (int)pow(STATES, 3);
  uint8_t rule[rule_size];
  uint32_t length[256] = {};
  size_t runs = 10;
  size_t size = 256;

  time_t t;
  srand((unsigned)time(&t));

  for (int c = 0; c < runs; ++c) {
    for (int n = 0; n < 256; n++) {
      for (int i = 0; i < rule_size; ++i) {
        rule[i] = (n & (1 << i)) ? 1: 0;
        /* printf("%"PRIu8, rule[i]); */
      }
      /* printf("\n"); */
      /* printf("%lu\n", rule_number(rule, RULE_SIZE)); */
      length[n] += write_to_file(size, rule_size, rule, 0, 1, STEPS);
    }
  }

  for (int n = 0; n < 256; n++) {
    printf("%i\t%f\n", n, ((float)length[n])/runs);
  }
  /* for (int rule = 0; rule < RULE_SIZE; rule++) { */
  /* printf("\rRule %i", rule); */
  /* write_to_file(rule, 0, 1); */
  /* } */
  return 0;
}

int main(int argc, char** argv)
{
  /* No arguments were passed */
  if (argc < 2) {
    printf("%s\n", help);
    return 0;
  }

  if (strcmp("2d", argv[1]) == 0) {
    return main_2d(argc, argv);
  }
  else if (strcmp("1d", argv[1]) == 0) {
    return main_1d(argc, argv);
  }
}
