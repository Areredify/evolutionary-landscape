#include <algorithm>
#include <iostream>
#include <numeric>
#include <iterator>
#include <fstream>
#include <ctime>
#include <thread>
#include <random>
#include <cmath>

#include "alg.h"

// function to output population info into a stream (stdio or file)
std::ostream &operator<<(std::ostream &stream, const population &pop)
{
    stream << pop.size << " " << pop.fitness_func.output_size << "\n";

    for (auto i = 0u; i < pop.size; i++)
    {
        for (auto j = 0u; j < pop.cur_pop[i].t_gene.size; j++)
            stream << pop.cur_pop[i].t_gene[j] << " ";
        stream << "\n";
    }

    return stream;
}

// async function to run a simulation
void simulation::run()
{
    const time_t t = clock();
    std::vector<time_t> times({t});

    // file with different statistics for the simulation (similarity, speed, ...)
    std::ofstream stats(out_path.append("statistics.txt"));

    // write initial populations to file
    std::ofstream(out_path.append("epoch0.txt")) << bacteria_pop << plant_pop;

    std::mt19937 random_generator((std::random_device())());
    unsigned epoch = 1;
    bacteria_pop.epoch = &epoch;
    plant_pop.epoch = &epoch;

    for (; epoch <= generations; epoch++)
    {
        evolve(random_generator);

        bacteria_pop.speed_scale += 1;
        plant_pop.speed_scale += 1;
        stats << epoch << " " << bacteria_pop.stats() << " " << plant_pop.stats() << " " << get_similarity() << " " << get_distance(random_generator) << std::endl;

        std::ostringstream oss;
        oss << "epoch" << epoch << ".txt";
        if ((epoch - 1) % 20 == 0)
            std::ofstream(std::filesystem::path(out_path).append(oss.str())) << bacteria_pop << plant_pop;

        times.push_back(clock());
        const auto dt = /*std::min(static_cast<size_t>(10),*/ times.size() - 1;
        const auto old_t = times[times.size() - dt - 1];
        std::cout << std::this_thread::get_id() << " " << epoch << " " << static_cast<double>(times.back() - old_t) / CLOCKS_PER_SEC / dt * (generations - epoch) << " " << static_cast<double>(times.back() - old_t) / CLOCKS_PER_SEC / dt << std::endl;
    }
}

double simulation::get_similarity() const
{
    auto [b_min, b_max] = bacteria_pop.min_max();
    auto [p_min, p_max] = plant_pop.min_max();

    const auto bp_min = itemwise_func(p_min, b_min, [](double x, double y) { return std::min(x, y); });
    const auto bp_max = itemwise_func(p_max, b_max, [](double x, double y) { return std::max(x, y); });

    auto b_discrete = bacteria_pop.labeling(bp_min, bp_max),
         p_discrete = plant_pop.labeling(bp_min, bp_max);

    double similarity = 0;
    for (auto j = 0u; j < b_discrete.size(); j++)
        similarity += sqrt(p_discrete[j] * b_discrete[j]);

    return sqrt(1 - similarity);
}

double simulation::get_distance(std::mt19937 &random_generator) const
{
    const unsigned samples = 50000;
    std::uniform_int_distribution<unsigned> b_gen(0, bacteria_pop.size - 1);
    std::uniform_int_distribution<unsigned> p_gen(0, plant_pop.size - 1);

    double res = 0;

    for (auto i = 0u; i < samples; i++)
        res += distance(bacteria_pop.cur_pop[b_gen(random_generator)].t_gene,
                        plant_pop.cur_pop[p_gen(random_generator)].t_gene);

    return res / samples;
}

double simulation::get_prob(const double mean)
{
    return 1 + 5 * pow(-cos(3) + cos(3 * std::min(mean, 1.0)), 3);
}

void simulation::evolve(std::mt19937 &random_generator)
{
    const unsigned bonus_reproduction = 20, bonus_bacteria = 10;
    std::vector<unsigned> bacteria_reproduction_p(bacteria_pop.size, 100);

    const auto pool_size = static_cast<unsigned>((bacteria_pop.size + plant_pop.size - 1) / plant_pop.size);

    std::vector<int> indices(bacteria_pop.size);
    for (auto i = 0u; i < indices.size(); i++)
        indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), random_generator);

    std::vector<double> distances(indices.size());

    for (auto i = 0u; i < distances.size(); i++)
        distances[i] = distance(bacteria_pop.cur_pop[indices[i]].t_gene,
                                plant_pop.cur_pop[i / pool_size].t_gene);

    // std::cout << *std::min_element(distances.begin(), distances.end()) << std::endl;

    std::vector<double> means;
    means.reserve(plant_pop.size);
    for (auto i = 0u; i < plant_pop.size; i++)
    {
        const auto l = i * pool_size, r = std::min(l + pool_size, static_cast<unsigned>(bacteria_pop.size));
        auto min_values = KMin(bonus_bacteria,
                               distances.begin() + l,
                               distances.begin() + r);
        double mean = 0;
        for (auto &j : min_values)
        {
            auto ind = std::distance(distances.begin(), j);
            ind = indices[ind];
            mean += *j;

            if (*j < 1)
                bacteria_reproduction_p[ind] += bonus_reproduction;
        }
        mean /= min_values.size();
        means.push_back(get_prob(mean));
    }

    std::discrete_distribution<> bd(bacteria_reproduction_p.begin(), bacteria_reproduction_p.end());

    const auto bacteria_prev = bacteria_pop.cur_pop;
    const auto plant_prev = plant_pop.cur_pop;

    bacteria_pop.advance_pop();
    plant_pop.advance_pop();

    for (auto i = 0u; i < bacteria_pop.size; i++)
        bacteria_pop.make_offspring(bacteria_prev + bd(random_generator), bacteria_pop.cur_pop + i, random_generator);

    std::discrete_distribution<> pd(means.begin(), means.end());
    for (auto i = 0u; i < plant_pop.size; i++)
        plant_pop.make_offspring(plant_prev + pd(random_generator), plant_pop.cur_pop + i, random_generator);
}
