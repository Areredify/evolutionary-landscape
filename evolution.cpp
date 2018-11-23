#include <algorithm>
#include <iostream>
#include <numeric>
#include <iterator>
#include <fstream>
#include <ctime>
#include <thread>

#include "alg.h"

std::ostream& operator<<(std::ostream& stream, const population& pop) {
    stream << pop.size << " " << pop.fitness_func.output_size << "\n";

    for (auto& i : pop.pop) {
        for (auto& j : i->t_gene)
            stream << j << " ";
        stream << "\n";
    }

    return stream;
}

void simulation::run() {
    const time_t t = clock();
    std::vector<time_t> times({ t });
    std::ofstream stats(std::experimental::filesystem::path(out_path).append("statistics.txt"));

    std::ofstream(std::experimental::filesystem::path(out_path).append("epoch0.txt")) << bacteria_pop << plant_pop;

    std::mt19937 random_generator((std::random_device())());
    unsigned epoch = 1;
    bacteria_pop.epoch = &epoch;
    plant_pop.epoch = &epoch;

    for (; epoch <= generations; epoch++) {
        evolve(random_generator);

        bacteria_pop.speed_scale += 1;
        plant_pop.speed_scale += 1;
        stats << epoch << " " << bacteria_pop.stats() << " " << plant_pop.stats() << " " << get_similarity() << " " << get_distance(random_generator) << std::endl;

        std::ostringstream oss;
        oss << "epoch" << epoch << ".txt";
        if ((epoch - 1) % 20 == 0)
            std::ofstream(std::experimental::filesystem::path(out_path).append(oss.str())) << bacteria_pop << plant_pop;


        times.push_back(clock());
        const auto dt = std::min(static_cast<size_t>(10), times.size() - 1);
        const size_t old_t = times[times.size() - dt - 1];
        std::cout << std::this_thread::get_id() << " " << epoch << " " << static_cast<double>(times.back() - old_t) / CLOCKS_PER_SEC / dt * (generations - epoch) << " " <<
                  static_cast<double>(times.back() - old_t) / CLOCKS_PER_SEC / dt << std::endl;
    }
}

double simulation::get_similarity() const {
    auto [b_min, b_max] = bacteria_pop.min_max();
    auto [p_min, p_max] = plant_pop.min_max();

    auto bp_min = itemwise_func(p_min, b_min, [](double x, double y) {return std::min(x, y); });
    auto bp_max = itemwise_func(p_max, b_max, [](double x, double y) {return std::max(x, y); });

    std::vector<double> b_discrete = bacteria_pop.labeling(bp_min, bp_max),
                        p_discrete = plant_pop.labeling(bp_min, bp_max);

    double similarity = 0;
    for (auto j = 0u; j < b_discrete.size(); j++)
        similarity += sqrt(p_discrete[j] * b_discrete[j]);

    return sqrt(1 - similarity); 
}

double simulation::get_distance(std::mt19937 & random_generator) const {
    const unsigned samples = 50000;
    std::uniform_int_distribution<unsigned> b_gen(0, bacteria_pop.size - 1);
    std::uniform_int_distribution<unsigned> p_gen(0, plant_pop.size - 1);

    double res = 0;

    for (auto i = 0u; i < samples; i++)
        res += distance(bacteria_pop.pop[b_gen(random_generator)]->t_gene,
                        plant_pop   .pop[p_gen(random_generator)]->t_gene);

    return res / samples;
}

double simulation::get_prob(double mean) {
    return 1 + 5 * pow(-cos(3) + cos(3 * std::min(mean, 1.0)), 3);
}

void simulation::evolve(std::mt19937 &random_generator) {
    const unsigned bonus_reproduction = 20, bonus_bacteria = 10;
    std::vector<unsigned> bacteria_reproduction_p(bacteria_pop.size, 100);

    const auto pool_size = static_cast<unsigned>((bacteria_pop.size + plant_pop.size - 1) / plant_pop.size);

    std::vector<int> indices(bacteria_pop.size);
    for (auto i = 0u; i < indices.size(); i++)
        indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), random_generator);

    std::vector<double> distances(indices.size());

    for (auto i = 0u; i < distances.size(); i++)
        distances[i] = distance(bacteria_pop.pop[indices[i]]->t_gene,
                                plant_pop.pop[i / pool_size]->t_gene);

    // std::cout << *std::min_element(distances.begin(), distances.end()) << std::endl;

    std::vector<double> means;
    means.reserve(plant_pop.size);
    for (auto i = 0u; i < plant_pop.size; i++) {
        const auto l = i * pool_size, r = std::min(l + pool_size, static_cast<unsigned>(bacteria_pop.size));
        auto min_values = KMin(bonus_bacteria, 
                               distances.begin() + l, 
                               distances.begin() + r);
        double mean = 0;
        for (auto &j : min_values) {
            auto ind = std::distance(distances.begin(), j);
            ind = indices[ind];
            mean += *j;

            if (*j < 1)
                bacteria_reproduction_p[ind] += bonus_reproduction;
        }
        mean /= min_values.size();
        means.push_back(get_prob(mean));
    }

    std::vector<std::shared_ptr<organism>> bacteria_res, plant_res;
    bacteria_res.reserve(bacteria_pop.size);
    plant_res.reserve(plant_pop.size);

    std::discrete_distribution<> bd(bacteria_reproduction_p.begin(), bacteria_reproduction_p.end());
    for (auto i = 0u; i < bacteria_pop.size; i++)
        bacteria_res.push_back(bacteria_pop.make_offspring(bacteria_pop.pop[bd(random_generator)], random_generator));

    std::discrete_distribution<> pd(means.begin(), means.end());
    for (auto i = 0u; i < plant_pop.size; i++)
        plant_res.push_back(plant_pop.make_offspring(plant_pop.pop[pd(random_generator)], random_generator));

    if (*bacteria_pop.epoch > bacteria_pop.ancestor_size)
        for (auto & i : bacteria_pop.pop)
           i->k_ancestor->ancestors.clear(), i->k_ancestor->k_ancestor = nullptr;

    if (*plant_pop.epoch > plant_pop.ancestor_size)
        for (auto & i : plant_pop.pop)
            i->k_ancestor->ancestors.clear(), i->k_ancestor->k_ancestor = nullptr;

    bacteria_pop.pop = std::move(bacteria_res);
    plant_pop.pop = std::move(plant_res);
}
