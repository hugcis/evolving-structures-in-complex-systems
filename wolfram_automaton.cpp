#include <iostream>
#include <sstream>
#include <fstream>
#include "zlib.h"

# define BITS 80
# define STEPS 800

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

unsigned char * printBits( int a[] , int bits_to_read )
{
  unsigned char * buf = (unsigned char *)malloc(bits_to_read * sizeof(unsigned char));
  for (int i = 0 ; i < bits_to_read ; i++) {
    buf[i] = (a[i/32] & (1 << (i%32))) != 0 ? '1': '0';
  }
  return buf;
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

int compress_memory(void *in_data, size_t in_data_size)
{
  const size_t BUFSIZE = 128 * 1024;
  uint8_t temp_buffer[BUFSIZE];

  z_stream strm;
  strm.zalloc = 0;
  strm.zfree = 0;
  strm.next_in = reinterpret_cast<uint8_t *>(in_data);
  strm.avail_in = in_data_size;
  strm.next_out = temp_buffer;
  strm.avail_out = BUFSIZE;

  deflateInit(&strm, Z_BEST_COMPRESSION);

  while (strm.avail_in != 0)
  {
    int res = deflate(&strm, Z_NO_FLUSH);
    assert(res == Z_OK);
    if (strm.avail_out == 0)
    {
      strm.next_out = temp_buffer;
      strm.avail_out = BUFSIZE;
    }
  }

  int deflate_res = Z_OK;
  while (deflate_res == Z_OK)
  {
    if (strm.avail_out == 0)
    {
      strm.next_out = temp_buffer;
      strm.avail_out = BUFSIZE;
    }
    deflate_res = deflate(&strm, Z_FINISH);
  }

  assert(deflate_res == Z_STREAM_END);
  int compressed_size = BUFSIZE - strm.avail_out;
  deflateEnd(&strm);

  return compressed_size;
}

void init_automat(int a[], int size)
{
  SetBit(a, BITS/2);
//   time_t t;
//   srand((unsigned) time(&t));
//   for (int i = 0 ; i < size ; i++) {
//     if (rand()%2 == 0) {
//       SetBit(a, i);
//     }
//   }
}

void write_to_file(int rule, int print_automaton)
{
  std::ofstream out_file;
  std::ostringstream fname;
  fname << "data/out" << rule << ".dat";
  out_file.open(fname.str());

  int A[BITS/32 + 1]  =  { };
  init_automat(A, BITS);
  for ( int i = 0 ; i < STEPS ; i++ ) {

    if (print_automaton == 1) { std::cout << printBits(A, BITS) << "\n"; }

    update_step(A, BITS, rule);
    if (i%30 == 0) {
      out_file << i << "    "
               << compress_memory(printBits(A, BITS), BITS) << "\n";
    }
  }
  out_file.close();
}

int main()
{
  for (int rule = 88 ; rule < 89 ; rule ++) {
    write_to_file(rule, 1);
  }
  return 0;
}
