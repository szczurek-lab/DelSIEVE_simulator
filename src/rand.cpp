//
// Created by senbaikang on 29.05.22.
//

#include "rand.h"

void RandomDevice::SetRandomDevice(unsigned long _seed) {
  this->seed = _seed;
  this->rd.seed(this->seed);
}

unsigned int GetRandDevSize() { return sizeof(RandomDevice); }

RandomDevice *AllocateMem(size_t len) {
  auto *tmp = (RandomDevice *)calloc(len, sizeof(RandomDevice));

  //  fprintf(stderr, "%p\n", &tmp[0]);
  //  fprintf(stderr, "%p\n", &tmp[1]);
  //  fprintf(stderr, "%p\n", &tmp[2]);
  //  fprintf(stderr, "%p\n", &tmp[3]);
  //  fprintf(stderr, "%p\n", &tmp[4]);

  return tmp;
}

void FreeMem(RandomDevice *rd) { free(rd); }

RandomDevice *GetRandDevElement(RandomDevice *head, size_t len) {
  return head + len;
}

void Call_C_SetRandomDevice(unsigned long seed, RandomDevice *rd) {
  rd->SetRandomDevice(seed);
}

unsigned long RandInteger(RandomDevice *rd) { return rd->rd(); }

void RandIntegers(int size, unsigned long *out, RandomDevice *rd) {
  for (size_t i = 0; i < size; i++) {
    out[i] = RandInteger(rd);
  }
}

double RandBeta(double mean, double var, RandomDevice *rd) {
  double sample_size, shape1, shape2, gamma1, gamma2, randBeta;

  /* assuming variance < mean (1-mean) */
  sample_size = (mean * (1.0 - mean) / var) - 1.0;

  shape1 = mean * sample_size;
  shape2 = (1.0 - mean) * sample_size;

  gamma1 = RandGamma(shape1, rd);
  gamma2 = RandGamma(shape2, rd);
  randBeta = gamma1 / (gamma1 + gamma2);

  return randBeta;
}

int RandBinomial(double prob, int numTrials, RandomDevice *rd) {
  std::binomial_distribution<int> dist(numTrials, prob);
  return dist(rd->rd);
}

double RandExponential(double lambda, RandomDevice *rd) {
  std::exponential_distribution<double> dist(lambda);
  return dist(rd->rd);
}

double RandGamma(double shape, RandomDevice *rd) {
  std::gamma_distribution<double> dist(shape);
  return dist(rd->rd);
}

int RandNegativeBinomial(double mean, double dispersion, RandomDevice *rd) {
  int poissonRand;
  double gammaRand;

  /* the RandomGamma function here has mean 1, so we need to scale it
   * ourselves
   */
  do {
    gammaRand = mean / dispersion * RandGamma(dispersion, rd);
  } while (fabs(gammaRand - 0) < 1e-16);
  poissonRand = RandPoisson(gammaRand, rd);

  return poissonRand;
}

int RandPoisson(double lambda, RandomDevice *rd) {
  std::poisson_distribution<int> dist(lambda);
  return dist(rd->rd);
}

double RandUniform(RandomDevice *rd) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rd->rd);
}

/* [from, to - 1] */
int RandUniformInt(int from, int to, RandomDevice *rd) {
  if (from >= to - 1)
  {
    fprintf(stderr, "Error! Unable to sample an integer between [%d,%d]", from, to - 1);
    exit(1);
  }

  std::uniform_int_distribution<int> dist(from, to - 1);
  return dist(rd->rd);
}

int RandBetaBinomial(double mean, double overdispersion, int numTrials,
                     RandomDevice *rd) {
  double alpha, beta, gamma1, gamma2, prob;
  int sum;

  alpha = mean * overdispersion;
  beta = overdispersion * (1 - mean);

  gamma1 = RandGamma(alpha, rd);
  gamma2 = RandGamma(beta, rd);
  prob = gamma1 / (gamma1 + gamma2);

  sum = 0;
  for (int i = 0; i < numTrials; i++) {
    if (RandUniform(rd) < prob)
      sum++;
  }
  return sum;
}

void Normalize(double *arr, size_t size) {
  double sum = 0.0;

  for (size_t i = 0; i < size; i++)
    sum += arr[i];

  if (sum == 1.0)
    return;

  for (size_t i = 0; i < size; i++)
    arr[i] = arr[i] / sum;
}

void RandDirichletMultinomial(double *means, double overdispersion,
                              int numTrials, int pivot, size_t size,
                              int *readCounts, RandomDevice *rd) {
  double gammas[size], probs[size], sum = 0.0;
  int badReads, readsLeft;

  Normalize(means, size);

  for (size_t i = 0; i < size; i++) {
    gammas[i] = RandGamma(means[i] * overdispersion, rd);
    sum += gammas[i];

    readCounts[i] = 0;
  }

  for (size_t i = 0; i < size; i++)
    probs[i] = gammas[i] / sum;

  readsLeft = numTrials;
  for (size_t i = 0; i < size; i++) {
    if (i != pivot) {
      badReads = RandBinomial(probs[i], numTrials, rd);

      if (readsLeft - badReads >= 0)
        readCounts[i] = badReads;
      else {
        badReads = readsLeft;
        readCounts[i] = badReads;
      }

      readsLeft -= badReads;

      if (readsLeft == 0)
        break;
    }
  }
  readCounts[pivot] = readsLeft;
}

double RandGaussian(double mean, double stddev, RandomDevice *rd) {
  std::normal_distribution<double> dist(mean, stddev);
  return dist(rd->rd);
}

void RandGaussianArr(size_t num, double mean, double stddev, RandomDevice *rd,
                     double *out, int forcePositive) {
  for (size_t i = 0; i < num; i++) {
    while (true) {
      out[i] = RandGaussian(mean, stddev, rd);

      if (forcePositive == 0 || out[i] > 0)
        break;
    }
  }
}
