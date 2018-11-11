#include <iostream>
#include <random>
#include <cmath>
#include <thread>
#include <cassert>

#include "alg.h"

double rnd0(std::mt19937 &random_generator) {
  static std::uniform_real_distribution<double> rnd(0.0, 1.0);

  return rnd(random_generator);
}


double rnd1(std::mt19937 &random_generator) {
  static std::uniform_real_distribution<double> rnd(-1.0, 1.0);

  return rnd(random_generator);
}

double nrm0(std::mt19937 &random_generator) {
  static std::normal_distribution<double> rnd(0, 1);

  return rnd(random_generator);
}

inline double sqr(const double x) {
  return x * x;
}

double distance(const std::vector<double> & v1,
                const std::vector<double> & v2)
{
    double sum = 0;
    //assert(v1.size() == v2.size());

    for (auto i = 0u; i < v1.size(); i++)
        sum += sqr(v1[i] - v2[i]);
    return sqrt(sum);
}


void print_gene(const std::vector<double> & v)
{
    for (auto & i : v)
        std::cout << i << " ";
    std::cout << "\n";
}

