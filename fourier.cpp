#define _USE_MATH_DEFINES
#include <cmath>

#include "alg.h"

unsigned fourier::index(const int o, const int i, const int f) const {
  return o * input_size * func_num + i * func_num + f;
}

fourier::fourier(const unsigned input_size, const unsigned output_size, const double power, const unsigned func_num)
    : fitness(input_size, output_size), func_num(func_num),
    transformer(input_size * output_size * func_num) {

    std::mt19937 random_generator((std::random_device())());

    for (auto i = 0u; i < output_size; i++)
        for (auto j = 0u; j < input_size; j++)
            for (auto k = 0u; k < func_num; k++)
                transformer[index(i, j, k)] = cell{ pow(k + 1, -power) * rnd1(random_generator),
                                                    2 * M_PI * k, 
                                                    2 * M_PI * rnd0(random_generator)};
}

std::vector<double> fourier::operator()(const std::vector<double> & gene) const {
    std::vector<double> res(output_size, 0);
    for (auto i = 0u; i < output_size; i++)
        for (auto j = 0u; j < input_size; j++)
            for (auto k = 0u; k < func_num; k++) {
                const cell & cur_cell = transformer[index(i, j, k)];
                res[i] += cur_cell.a * cos(cur_cell.k * gene[j] + cur_cell.b);
            }
    return res;
}