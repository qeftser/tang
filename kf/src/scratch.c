
#include "numerical_basics.h"
#include <stdio.h>

DEFMAT(10,10);

int main(void) {

   mat10x10 S = MAT(10,10);
   mat10x10 Si = MAT(10,10);
   mat10x10 SS = MAT(10,10);

   for (int i = 0; i < 100; ++i)
      ((matbase *)&S)->vals[i] = (double)rand() / RAND_MAX;


   matshw(&S);
   printf("%f\n",matdet(&S));

   matshw(matmul(matinv(&S,&Si),&S,&SS));

   return 0;
}
