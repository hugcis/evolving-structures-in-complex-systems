/**
 * Program that generates ppm file from the step file.
 * Compile with
 * gcc step_to_ppm.c -o step_to_ppm -lnetpbm
 */

#include <stdio.h>
#include <string.h>
#include <netpbm/pam.h>

#define CHANNELS 3

const int palette_rgb[23][CHANNELS] = {
                                {255, 255, 255},
                                {0, 0, 0},
                                {255, 0, 0},     /* Red */
                                {0, 255, 0},     /* Green */
                                {0, 0, 255},     /* Blue */
                                {127, 0, 255},   /* Purple */
                                {255, 69, 0},    /* Red-Orange */
                                {255, 255, 0},   /* Yellow */
                                {13, 152, 186},  /* Blue-green */
                                {255, 0, 0},
                                {204, 255, 0},
                                {0, 255, 101},
                                {0, 101, 255},
                                {203, 0, 255},
                                {255, 100, 51},
                                {42, 75, 215},
                                {29, 105, 20},
                                {129, 74, 25},
                                {129, 38, 192},
                                {41, 208, 208},
                                {255, 146, 51},
                                {233, 222, 187},
                                {255, 205, 243}};

int main(int argc, const char** argv)
{
  struct pam outpam;
  int size = atoi(argv[2]);
  int states = atoi(argv[3]);
  int row, c;
  FILE* input = fopen(argv[1], "r");
  uint8_t automaton[size * size];

  int i = 0;
  while ((c = fgetc(input)) != EOF) {
    if (c != '\n') {
      automaton[i++] = c - '0';
    }
  }

  pm_proginit(&argc, argv);

  outpam.size = 10000;
  outpam.len = 10000;
  outpam.maxval = 255;
  strcpy(outpam.tuple_type, PAM_PPM_TUPLETYPE);
  outpam.format = PAM_FORMAT;
  outpam.height = size;
  outpam.width = size;
  outpam.file = stdout;
  outpam.depth = CHANNELS;

  tuple** arr = pnm_allocpamarray(&outpam);
  int color_offset;
  if (states == 2) {
    color_offset = 0;
  } else if (states == 3) {
    color_offset = 2;
  } else if (states == 4) {
    color_offset = 5;
  } else {
    color_offset = 9;
  }

  for (row = 0; row < outpam.height; ++row) {
    int column;
    for (column = 0; column < outpam.width; ++column) {
      unsigned int depth;
      for (depth = 0; depth < CHANNELS; ++depth) {
        arr[row][column][depth] =
          palette_rgb[automaton[row * size + column] +
                      color_offset][depth];
      }
    }
  }
  pnm_writepam(&outpam, arr);
  pnm_freepamarray(arr, &outpam);
}
