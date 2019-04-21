#include <iostream>
#include <random>
#include <cmath>
#include <thread>
#include <cassert>

#include "alg.h"

// gives random number uniformly distributed from 0 to 1
double rnd0(std::mt19937 &random_generator) {
    static std::uniform_real_distribution<double> rnd(0.0, 1.0);
    
    return rnd(random_generator);
}


// gives random number uniformly distributed from -1 to 1
double rnd1(std::mt19937 &random_generator) {
    static std::uniform_real_distribution<double> rnd(-1.0, 1.0);
    
    return rnd(random_generator);
}

// gives random number normally distributed with mean 0 and dispersion 1
double nrm0(std::mt19937 &random_generator) {
    static std::normal_distribution<double> rnd(0, 1);
    
    return rnd(random_generator);
}

// random gene with coordinates each distributed uniformly from -1 to 1
gene_c rnd1_vec(unsigned len, std::mt19937 &random_generator) {
    gene_c res;
    res.size = len;

    for (auto i = 0u; i < len; i++)
        res[i] = rnd1(random_generator);

    return res;
}

double sqr(const double x) {
    return x * x;
}


// euclidian distance between two vectors
double distance(const gene_c & v1,
                const gene_c & v2) {
    double sum = 0;
    //assert(v1.size() == v2.size());

    for (auto i = 0u; i < v1.size; i++)
        sum += sqr(v1[i] - v2[i]);
    return sqrt(sum);
}


// squared euclidian distance between two vectors
double distance2(const gene_c & v1,
                 const gene_c & v2) {
    double sum = 0;
    //assert(v1.size() == v2.size());

    for (auto i = 0u; i < v1.size; i++)
        sum += sqr(v1[i] - v2[i]);
    return sum;
}