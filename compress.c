#include <stdlib.h>
#include <assert.h>
#include "compress.h"
#include "zlib.h"

int compress_memory_size(void *in_data, size_t in_data_size)
{
  const size_t BUFSIZE = 128 * 1024;
  uint8_t temp_buffer[BUFSIZE];

  z_stream strm;
  strm.zalloc = 0;
  strm.zfree = 0;
  strm.next_in = (uint8_t *)(in_data);
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
