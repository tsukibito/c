#ifndef MMUTIL
#define MMUTIL

#include "matrix_vector.h"
#define MM_MAX_LINE_BUF 256 

matrix mm_matrix(char *fileName, int *n, int *m, int start);
void mm_text(char *inputName, char *outputName);

#endif
