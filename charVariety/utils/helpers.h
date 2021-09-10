// charVariety/utils/helpers.h

#ifndef CHARVARIETY_HELPERS
#define CHARVARIETY_HELPERS

#include "dataStructure.h"

Pair DTX(Pair pr);
Pair DTY(Pair pr);
Pair DTX1(Pair pr);
Pair DTY1(Pair pr);
Pair action(Pair aa, int* maps, int nIter);

double avgLogExpansion(Matrix* A, int nMatrix, double angle);
double minAvgExpansionNewton(Matrix* A, int nMatrix);

#endif