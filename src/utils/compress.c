#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "compress.h"
#include "zlib.h"

void compress_in_memory(void* in_data, size_t in_data_size,
                        size_t bufsize, uint8_t* temp_buffer,
                        z_stream* strm)
{
  strm->zalloc = 0;
  strm->zfree = 0;
  strm->next_in = (uint8_t *)(in_data);
  strm->avail_in = in_data_size;
  strm->next_out = temp_buffer;
  strm->avail_out = bufsize;
  strm->data_type = Z_TEXT;

  deflateInit(strm, Z_BEST_COMPRESSION);

  while (strm->avail_in != 0)
    {
      int res = deflate(strm, Z_NO_FLUSH);
      assert(res == Z_OK);
      if (strm->avail_out == 0)
        {
          strm->next_out = temp_buffer;
          strm->avail_out = bufsize;
        }
    }

  int deflate_res = Z_OK;
  while (deflate_res == Z_OK)
    {
      if (strm->avail_out == 0)
        {
          strm->next_out = temp_buffer;
          strm->avail_out = bufsize;
        }
      deflate_res = deflate(strm, Z_FINISH);
    }

  assert(deflate_res == Z_STREAM_END);
}


int compress_memory_size(void *in_data, size_t in_data_size)
{
  const size_t BUFSIZE = 256 * 256 * 256;
  uint8_t* temp_buffer = malloc(sizeof(char) * BUFSIZE);
  z_stream strm = {};
  compress_in_memory(in_data, in_data_size, BUFSIZE, temp_buffer, &strm);
  int compressed_size = BUFSIZE - strm.avail_out;
  deflateEnd(&strm);
  free(temp_buffer);

  return compressed_size;
}

void compress_rule(char* rule_buf, uint8_t* out_buf, size_t buf_size)
{
  z_stream strm = {};
  compress_in_memory(rule_buf, strlen(rule_buf), buf_size, out_buf, &strm);
}
