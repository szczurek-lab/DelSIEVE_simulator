#include <climits>
//
// Created by Senbai Kang on 7/23/20.
//

#include "statisticalTools.h"
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <dlib/global_optimization.h>
#include <vector>

/* declare and define data structures */
struct OptimizeBetaBinMeanOverDis {
  std::vector<std::pair<unsigned, unsigned>> const &counts;

  explicit OptimizeBetaBinMeanOverDis(
      std::vector<std::pair<unsigned, unsigned>> const &counts_)
      : counts(counts_){};

  double operator()(const dlib::matrix<double, 0, 1> &x) const {
    double result = 0;
    for (size_t cell = 0; cell < this->counts.size(); ++cell) {
      result += logBetaBinPDF(this->counts[cell].first,
                              this->counts[cell].second, x(0), x(1));
    }
    return result;
  };
};

struct OptimizeBetaBinMeanOverDisDerivates {
  std::vector<std::pair<unsigned, unsigned>> const &counts;

  explicit OptimizeBetaBinMeanOverDisDerivates(
      std::vector<std::pair<unsigned, unsigned>> const &counts_)
      : counts(counts_) {}

  dlib::matrix<double> operator()(const dlib::matrix<double, 0, 1> &x) const {
    double mean = x(0);
    double overDis = x(1);
    dlib::matrix<double, 0, 1> res = {0, 0};

    double temp = 0;
    unsigned counter = 0;
    for (size_t cell = 0; cell < this->counts.size(); ++cell) {
      unsigned k = this->counts[cell].first;
      unsigned n = this->counts[cell].second;
      temp += overDis * boost::math::digamma(k + mean * overDis) -
              overDis * boost::math::digamma(n - k + overDis - overDis * mean);
      ++counter;
    }
    res(0) =
        counter * (-overDis * boost::math::digamma(mean * overDis) +
                   overDis * boost::math::digamma(overDis - overDis * mean)) +
        temp;

    temp = 0;
    for (size_t cell = 0; cell < this->counts.size(); ++cell) {
      unsigned k = this->counts[cell].first;
      unsigned n = this->counts[cell].second;
      temp += mean * boost::math::digamma(k + mean * overDis) +
              (1.0 - mean) *
                  boost::math::digamma(n - k + overDis - overDis * mean) -
              boost::math::digamma(n + overDis);
    }
    res(1) = counter * (boost::math::digamma(overDis) -
                        mean * boost::math::digamma(mean * overDis) -
                        (1.0 - mean) *
                            boost::math::digamma(overDis - overDis * mean)) +
             temp;

    return res;
  }
};

struct OptimizeBetaBinOverDis {
  std::vector<std::pair<unsigned, unsigned>> const &counts;
  double meanFilter;

  OptimizeBetaBinOverDis(
      std::vector<std::pair<unsigned, unsigned>> const &counts_,
      double meanFilter_)
      : counts(counts_), meanFilter(meanFilter_){};

  double operator()(const dlib::matrix<double, 0, 1> &x) const {
    double result = 0;
    for (size_t cell = 0; cell < this->counts.size(); ++cell) {
      result += logBetaBinPDF(this->counts[cell].first,
                              this->counts[cell].second, meanFilter, x(0));
    }
    return result;
  };
};

struct OptimizeBetaBinOverDisDerivates {
  std::vector<std::pair<unsigned, unsigned>> const &counts;
  double meanFilter;

  OptimizeBetaBinOverDisDerivates(
      std::vector<std::pair<unsigned, unsigned>> const &counts_,
      double meanFilter_)
      : counts(counts_), meanFilter(meanFilter_){};

  dlib::matrix<double> operator()(const dlib::matrix<double, 0, 1> &x) const {
    double mean = this->meanFilter;
    double overDis = x(0);
    dlib::matrix<double, 0, 1> res = {0};

    double temp = 0;
    unsigned counter = 0;
    for (size_t cell = 0; cell < this->counts.size(); ++cell) {
      unsigned k = this->counts[cell].first;
      unsigned n = this->counts[cell].second;
      temp += mean * boost::math::digamma(k + mean * overDis) +
              (1.0 - mean) *
                  boost::math::digamma(n - k + overDis - overDis * mean) -
              boost::math::digamma(n + overDis);
      ++counter;
    }
    res(0) = counter * (boost::math::digamma(overDis) -
                        mean * boost::math::digamma(mean * overDis) -
                        (1.0 - mean) *
                            boost::math::digamma(overDis - overDis * mean)) +
             temp;

    return res;
  }
};

/* implement methods */
double logBetaBinCountsTerm(double sup, double cov) {
  return std::lgamma(cov + 1.0) - std::lgamma(sup + 1.0) -
         std::lgamma(cov - sup + 1.0);
}

double logBetaBinMixedTerm(double sup, double cov, double mean,
                           double overDis) {
  return std::lgamma(sup + mean * overDis) +
         std::lgamma(cov - sup + overDis * (1.0 - mean)) -
         std::lgamma(cov + overDis);
}

double logBetaBinParamsTerm(double mean, double overDis) {
  return std::lgamma(overDis) - std::lgamma(mean * overDis) -
         std::lgamma(overDis * (1.0 - mean));
}

double logBetaBinPDF(double sup, double cov, double mean, double overDis) {
  if (cov == 0)
    return 0;

  return logBetaBinCountsTerm(sup, cov) +
         logBetaBinMixedTerm(sup, cov, mean, overDis) +
         logBetaBinParamsTerm(mean, overDis);
}

int likelihoodRatioTest(CellStrExtd const *cellStr, int cellNum, int healthyTip,
                        size_t site, int altNuc, int cellGroupMode,
                        double meanFilter) {
  std::vector<std::pair<unsigned, unsigned>> testCounts;
  for (size_t cell = 0; cell != cellNum; cell++) {
    if (cell != healthyTip &&
        ((cellGroupMode == 0 &&
          cellStr[cell].site[site].readCount[altNuc] > 0) ||
         cellGroupMode == 1)) {
      testCounts.emplace_back(cellStr[cell].site[site].readCount[altNuc],
                              cellStr[cell].site[site].numReads);
    }
  }

  /* fit the relaxed model */
  dlib::matrix<double, 0, 1> startingPointMeanOverDis = {0.05, 5.0};
  dlib::matrix<double, 0, 1> dLibMinMeanOverDis = {0.001, 0.1};
  dlib::matrix<double, 0, 1> dLibMaxMeanOverDis = {0.999, 10000.0};
  OptimizeBetaBinMeanOverDis optBetaBinMeanOverDis(testCounts);
  OptimizeBetaBinMeanOverDisDerivates optBetaBinDerMeanOverDis(testCounts);
  double resultMeanOverDis = dlib::find_max_box_constrained(
      dlib::bfgs_search_strategy(), // Use BFGS search algorithm
      dlib::objective_delta_stop_strategy(
          1e-7), // Stop when the change in rosen() is less than 1e-7
      optBetaBinMeanOverDis, optBetaBinDerMeanOverDis, startingPointMeanOverDis,
      dLibMinMeanOverDis, dLibMaxMeanOverDis);

  /* fit the constrained model */
  dlib::matrix<double, 0, 1> startingPointOverDis = {2.0};
  dlib::matrix<double, 0, 1> dLibMinOverDis = {0.1};
  dlib::matrix<double, 0, 1> dLibMaxOverDis = {10000.0};
  OptimizeBetaBinOverDis optBetaBinOverDis(testCounts, meanFilter);
  OptimizeBetaBinOverDisDerivates optBetaBinDerOverDis(testCounts, meanFilter);
  double resultOverDis = dlib::find_max_box_constrained(
      dlib::bfgs_search_strategy(), // Use BFGS search algorithm
      dlib::objective_delta_stop_strategy(
          1e-7), // Stop when the change in rosen() is less than 1e-7
      optBetaBinOverDis, optBetaBinDerOverDis, startingPointOverDis,
      dLibMinOverDis, dLibMaxOverDis);

  /* get p-value */
  double pValue;
  if (resultOverDis > resultMeanOverDis) // if true the results are within the
                                         // optimization limit
  {
    pValue = 1;
  } else {
    double testStat = -2 * (resultOverDis - resultMeanOverDis);
    boost::math::chi_squared mydist(1);
    pValue = 1 - boost::math::cdf(mydist, testStat);
  }

  if (pValue > 0.05 || startingPointMeanOverDis(0) >= meanFilter) {
    /* SNV */
    return 1;
  } else {
    /* not SNV */
    return 0;
  }
}