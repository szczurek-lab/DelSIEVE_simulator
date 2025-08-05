//
// Created by senbaikang on 29.05.22.
//

#ifndef MYSIMULATOR_RAND_H
#define MYSIMULATOR_RAND_H

#ifdef __cplusplus
#include <cstdlib>
#include <random>
#else
#include <stdlib.h>
#endif

typedef struct RandomDevice {
  unsigned long seed;

#ifdef __cplusplus
  std::mt19937_64 rd;

  void SetRandomDevice(unsigned long _seed);
#endif

} RandomDevice;

#ifdef __cplusplus
extern "C" {
#endif

unsigned int GetRandDevSize();
RandomDevice *AllocateMem(size_t len);
void FreeMem(RandomDevice *rd);
RandomDevice *GetRandDevElement(RandomDevice *head, size_t len);

void Call_C_SetRandomDevice(unsigned long seed, RandomDevice *rd);
unsigned long RandInteger(RandomDevice *rd);
void RandIntegers(int size, unsigned long *out, RandomDevice *rd);
double RandBeta(double mean, double var, RandomDevice *rd);
int RandBinomial(double prob, int numTrials, RandomDevice *rd);
double RandExponential(double lambda, RandomDevice *rd);
double RandGamma(double shape, RandomDevice *rd);
int RandNegativeBinomial(double mean, double dispersion, RandomDevice *rd);
int RandPoisson(double lambda, RandomDevice *rd);
double RandUniform(RandomDevice *rd);
int RandUniformInt(int from, int to, RandomDevice *rd);

int RandBetaBinomial(double mean, double overdispersion, int numTrials,
                     RandomDevice *rd);
void Normalize(double *arr, size_t size);
void RandDirichletMultinomial(double *means, double overdispersion,
                              int numTrials, int pivot, size_t size,
                              int *readCounts, RandomDevice *rd);
double RandGaussian(double mean, double stddev, RandomDevice *rd);
void RandGaussianArr(size_t num, double mean, double stddev, RandomDevice *rd,
                     double *out, int forcePositive);

#ifdef __cplusplus
}
#endif

#endif // MYSIMULATOR_RAND_H
