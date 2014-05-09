#include <math.h>
#include <stdio.h>
#include <time.h>

void diagMul(int size, float* xx, float* res) {
  /* clock_t cstart = clock(); */
  int dim = sqrt(size);
  int _xRes = sqrt(size);
  int y;
  int x;
  for (y = 0; y < dim; y++) {
    for (x = 0; x < dim; x++) {
      if(x>0) {
        res[x+y*_xRes] += - xx[(x-1)+y*_xRes];
      }
      if(x < _xRes-1) {
        res[x+y*_xRes] += - xx[(x+1)+y*_xRes];
      }
      if(y>0) {
        res[x+y*_xRes] += - xx[x+(y-1)*_xRes];
      }
      if(y < _xRes-1) {
        res[x+y*_xRes] += - xx[x+(y+1)*_xRes];
      }
    }
  }
  /* clock_t cend = clock(); */
  /* printf("a %f\n", (double)cend - (double)cstart); */
  /* printf("a %i\n", size); */
}

void jacobi(int size, float* x, float* b, float* res) {
  /* clock_t cstart = clock(); */
  float w = 2/3.0;
  float Rx[size];
  diagMul(size, x, Rx);
  int i;
  for(i = 0; i < size; i++) {
    res[i] = w / 4.0 * (b[i] - Rx[i]) + (1 - w) * x[i];
  }
  /* clock_t cend = clock(); */
  /* printf("b %f\n", (double)cend - (double)cstart); */
  /* printf("b %i\n", size); */
}
