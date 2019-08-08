#include <stdlib.h>

#ifndef COMPRESS_H   /* Include guard */
#define COMPRESS_H

int compress_memory_size(void*, size_t);

void compress_rule(char* rule_buf, uint8_t* out_buf, size_t buf_size);

#endif // COMPRESS_H
