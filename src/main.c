#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "2d_automaton.h"

int main(int argc, char** argv)
{
  char rule_buf[GRULE_SIZE + 1];
  uint8_t rule_array[GRULE_SIZE];

  /* Write steps for a given rule  */
  if (argc == 2) {
    build_rule_from_args(rule_array, rule_buf, argv);
    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("Rule %s\n", rule_buf);

    process_rule(rule_array, rule_buf, 1, 1, 0);
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

    process_rule(rule_array, rule_buf, 0, 1, 0);
  }
  printf("\n");
  return 0;
}
