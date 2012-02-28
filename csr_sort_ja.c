#include <stdlib.h>

static int int_cmp(const void *a, const void *b)
{
  return *((int *)a) - *((int *)b);
}


void csr_sort_ja(int m, int *ia, int *ja)
{
  int     i;

  for (i=0; i<m; i++) {
    qsort(ja + ia[i] - 1, ia[i+1] - ia[i], sizeof(int), &int_cmp);
  }
}
