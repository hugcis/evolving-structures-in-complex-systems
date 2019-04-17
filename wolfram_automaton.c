#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "compress.h"

# define BITS 80
# define STEPS 80

void  SetBit( int A[],  int k )
{
  A[k/32] |= 1 << (k%32);  // Set the bit at the k-th position in A[i]
}

void  ClearBit( int A[],  int k )
{
  A[k/32] &= ~(1 << (k%32));
}

int TestBit( int A[],  int k )
{
  return ( (A[k/32] & (1 << (k%32) )) != 0 ) ;
}

unsigned char * print_bits_spaced( int a[] , int bits_to_read, int spaced )
{
  int n_bits = spaced == 1 ? 2*bits_to_read - 1: bits_to_read;
  int buf_index;

  unsigned char * buf = (unsigned char *)malloc(n_bits * sizeof(uint8_t));
  for (int i = 0 ; i < bits_to_read ; i++) {
    buf_index = spaced == 1 ? 2*i : i;
    buf[buf_index] = (a[i/32] & (1 << (i%32))) != 0 ? '1': '0';
    if (spaced == 1 & buf_index + 1 < n_bits) { buf[buf_index + 1] = ' '; }
  }
  return buf;
}

unsigned char * printBits( int a[], int bits_to_read )
{
  return print_bits_spaced(a, bits_to_read, 0);
}


void update_step(int base[], int len, int rule)
{
  int * A = (int *)malloc(len * sizeof(int));
  memcpy(A, base, len * sizeof(int));

  for (int r = 0 ; r < len ; r ++ ) {
    int c = 0;
    if ( TestBit(base, (r - 1 + len)%len) ) { c |= (1 << 2); }
    if ( TestBit(base, r) ) { c |= (1 << 1); }
    if ( TestBit(base, (r + 1)%len) ) { c |= (1 << 0); }
    if ( rule & (1 << c) ) {
      SetBit(A, r);
    } else {
      ClearBit(A, r);
    }
  }

  memcpy(base, A, len * sizeof(int));
}

void init_automat(int a[], int size)
{
  /* SetBit(a, BITS/2); */
  time_t t;
  srand((unsigned) time(&t));
  /* for (int i = 0 ; i < size ; i++) { */
    /* if (rand()%2 == 0) { */
      /* SetBit(a, i); */
    /* } */
  /* } */
  SetBit(a, 0);
  for (int i = 1 ; i < size/5 ; i++) {
    if (rand()%2 == 0) {
      SetBit(a, i);
    }
  }
}

void write_to_file(int rule, int print_automaton, int write)
{
  FILE *out_file;
  char * fname;
  asprintf(&fname, "data/out%i.dat", rule);
  out_file = fopen(fname, "w+");

  FILE *out_steps_file;
  char * out_steps_fname;
  asprintf(&out_steps_fname, "steps/out%i.steps", rule);
  out_steps_file = fopen(out_steps_fname, "w+");

  int A[BITS/32 + 1]  =  { };
  init_automat(A, BITS);

  for ( int i = 0 ; i < STEPS ; i++ ) {

    if (i%50 == 0 & write == 1) {
      FILE *steps_file;
      char * steps_fname;
      asprintf(&steps_fname, "steps/out%i_%i.step", rule, i);
      steps_file = fopen(steps_fname, "w+");
      fputs((char *)printBits(A, BITS), steps_file);
      fclose(steps_file);
    }

    fprintf(out_steps_file, "%s\n", print_bits_spaced(A, BITS, 1));
    if (print_automaton == 1) { printf("%s\n", printBits(A, BITS)); }

    update_step(A, BITS, rule);
    if (i%30 == 0) {
      fprintf(out_file, "%i    %i\n", i,
              compress_memory_size(printBits(A, BITS), BITS));
    }
  }
  fclose(out_file);
  fclose(out_steps_file);
}

int main()
{
  for (int rule = 0 ; rule < 256 ; rule ++) {
    write_to_file(rule, 0, 0);
  }
  return 0;
}
