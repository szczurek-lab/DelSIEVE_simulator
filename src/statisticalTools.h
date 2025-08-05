//
// Created by Senbai Kang on 7/23/20.
//

#ifndef MYSIMULATOR_STATISTICALTOOLS_H
#define MYSIMULATOR_STATISTICALTOOLS_H

#include "extraDefinitions.h"

#ifdef __cplusplus
#include <cstdlib>
extern "C" {
#else
#include <stdlib.h>
#endif

double logBetaBinCountsTerm(double, double);
double logBetaBinMixedTerm(double, double, double, double);
double logBetaBinParamsTerm(double, double);
double logBetaBinPDF(double, double, double, double);

int likelihoodRatioTest(CellStrExtd const *, int, int, size_t, int, int,
                        double);

#ifdef __cplusplus
}
#endif

#endif // MYSIMULATOR_STATISTICALTOOLS_H
