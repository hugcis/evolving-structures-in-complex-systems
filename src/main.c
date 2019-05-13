#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "2d_automaton.h"
#include "wolfram_automaton.h"

#define STEPS 1000

const char help[] = "Use with either 2d or 1d as first argument";

int main_2d(int argc, char** argv)
{
  const uint64_t grule_size = (int) pow(STATES, NEIGH_SIZE + 1);
  char rule_buf[grule_size + 1];
  uint8_t rule_array[grule_size];

  /* Write steps for a given rule  */
  if (argc == 3) {
    build_rule_from_args(grule_size, rule_array, rule_buf, argv[2]);
    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("Rule %s\n", rule_buf);

    process_rule(grule_size, rule_array, rule_buf, 1, 1, 0, STEPS);
    return 0;
  }


  /* Generate compression plots for many rules */
  time_t t;
  srand((unsigned)time(&t));
  FILE* dic_file;
  char* fname;

  for (int i = 0; i < 1000; ++i) {
    /* generate_totalistic_rule(rule_array, rule_buf); */
    generate_general_rule(grule_size, rule_array, rule_buf);

    asprintf(&fname, "data_2d_%i/%lu.map", STATES, hash(rule_buf));
    dic_file = fopen(fname, "w+");
    fprintf(dic_file, "%s", rule_buf);
    fclose(dic_file);

    sprintf(rule_buf, "%lu", hash(rule_buf));
    printf("%i: Rule %s\n", i, rule_buf);
    fflush(stdout);

    process_rule(grule_size, rule_array, rule_buf, 0, 1, 0, STEPS);
  }
  printf("\n");
  return 0;
}

int main_1d(int argc, char** argv)
{
  extern char *optarg;
  extern int optind;
  int c, err = 0;
  size_t size = 256, neighbor = 3;
  int states = 2;
  struct Options1D options;
  options.init = ONE;
  options.timesteps = 256;
  options.write = NO_WRITE;

  char usage[] = "%s 1d [-s size] [-t timesteps]"
                       "[-n n_states] [-w neighborhood width]"
                       "[-i (one|random|random_small)] [-o]";
  char invalid_init[] = "Invalid value \"%s\" for init option."
                        " Must be one of \"one\", \"random\","
                        " \"random_small\"\n";

  while ((c = getopt(argc - 1, &argv[1], "s:t:n:w:i:o")) != -1)
    switch (c) {
    case 's':
      size = atoi(optarg);
      break;
    case 't':
      options.timesteps = atoi(optarg);
      break;
    case 'n':
      states = atoi(optarg);
      break;
    case 'w':
      neighbor = atoi(optarg);
      break;
    case 'o':
      options.write = WRITE_STEP;
      break;
    case 'i':
      if (strcmp("one", optarg) == 0) {
        options.init = ONE;
      }
      else if (strcmp("random", optarg) == 0) {
        options.init = RANDOM;
      }
      else if (strcmp("random_small", optarg) == 0) {
        options.init = RAND_SMALL;
      }
      else {
        printf(invalid_init, optarg);
        exit(1);
      }
      break;
    case '?':
      err = 1;
      break;
    }

  if (err) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  size_t rule_size = (int)pow(states, neighbor);
  uint8_t rule[rule_size];

  time_t t;
  srand((unsigned)time(&t));

  for (int n = 0; n < (int)pow(states, rule_size); n++) {
    for (int i = 0; i < rule_size; ++i) {
      rule[i] = (n & (1 << i)) ? 1: 0;
    }
    printf("%lu\t", rule_number(rule_size, rule));
    write_to_file(size, rule_size, rule, 0, 1, &options);
    printf("\n");
  }

  return 0;
}

int main(int argc, char** argv)
{
  char base_usage[] = "%s [\"1d\" or \"2d\"]";
  /* No first argument was passed */
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
  else {
    printf("%s\n", help);
    fprintf(stderr, base_usage, argv[0]);
  }
}
