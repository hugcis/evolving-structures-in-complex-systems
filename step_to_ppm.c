/* Program that generates ppm file from the step file.
   Compile with
   gcc step_to_ppm.c -o step_to_ppm -lnetpbm
 */

#include <stdio.h>
#include <string.h>
#include <netpbm/pam.h>

const int palette_rgb[14][3] = {
                                {255, 255, 255},
                                {0, 0, 0},
                                {255, 0, 0},
                                {0, 255, 0},
                                {0, 0, 255},
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
  outpam.depth = 3;

  tuple** arr = pnm_allocpamarray(&outpam);
  int color_offset = states == 2 ? 0 : 2;

  for (row = 0; row < outpam.height; ++row) {
    int column;
    for (column = 0; column < outpam.width; ++column) {
      unsigned int depth;
      for (depth = 0; depth < outpam.depth; ++depth) {
        arr[row][column][depth] =
          palette_rgb[automaton[row * size + column] +
                      color_offset][depth];
      }
    }
  }
  pnm_writepam(&outpam, arr);
  pnm_freepamarray(arr, &outpam);
}
