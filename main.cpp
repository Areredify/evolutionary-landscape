#include <iostream>
#include <cstdlib>

#include "ctpl_stl.h"

#include "alg.h"

std::pair<org_vec, org_vec> random_haplotypes(const simulation_params& sp,
                                             const fourier& bacteria_fitness,
                                             const fourier& plant_fitness,
                                             const unsigned ancestor_size,
                                             const double mn, const double mx) {
    org_vec b_pop, p_pop;
    std::mt19937 random_generator(std::random_device{}());

    while (true) {
        std::vector<double> b_haplotype(sp.b_gene_size),
                            p_haplotype(sp.p_gene_size);

        for (auto & i : b_haplotype)
            i = rnd1(random_generator);

        for (auto & i : p_haplotype)
            i = rnd1(random_generator);

        const auto d = distance(bacteria_fitness(b_haplotype), plant_fitness(p_haplotype));

        if (d > mn && d < mx) {
            auto bacteria = std::make_shared<organism>(b_haplotype);
            bacteria->ancestors.resize(1 + static_cast<unsigned>(log2(ancestor_size) + 0.01));
            bacteria->t_gene = bacteria_fitness(b_haplotype);
        
            auto plant = std::make_shared<organism>(p_haplotype);
            plant->ancestors.resize(1 + static_cast<unsigned>(log2(ancestor_size) + 0.01));
            plant->t_gene = plant_fitness(p_haplotype);
                                                                                                                                                
            for (auto i = 0u; i < sp.b_pop_size; i++) {
                b_pop.push_back(bacteria);
            }
            for (auto i = 0u; i < sp.p_pop_size; i++) {
                p_pop.push_back(plant);
            }
            return {b_pop, p_pop};
        }
    }
}

std::pair<org_vec, org_vec> full_random(const simulation_params& sp,
                                        const fourier& bacteria_fitness,
                                        const fourier& plant_fitness,
                                        const unsigned ancestor_size) {
    org_vec b_pop, p_pop;
    std::mt19937 random_generator(std::random_device{}());
    for (auto i = 0u; i < sp.b_pop_size; i++) {
        auto bacteria = std::make_shared<organism>(rnd1_vec(bacteria_fitness.input_size, random_generator));
        bacteria->ancestors.resize(1 + static_cast<unsigned>(log2(ancestor_size) + 0.01));
        bacteria->t_gene = bacteria_fitness(bacteria->gene);
        b_pop.push_back(bacteria);
    }
    for (auto i = 0u; i < sp.p_pop_size; i++) {
        auto plant = std::make_shared<organism>(rnd1_vec(plant_fitness.input_size, random_generator));
        plant->ancestors.resize(1 + static_cast<unsigned>(log2(ancestor_size) + 0.01));
        plant->t_gene = plant_fitness(plant->gene);
        p_pop.push_back(plant);
    }

    return {b_pop, p_pop};
}

void run_simulation(int id, const simulation_params& sp) {
    const fourier bacteria_fitness(sp.b_gene_size, sp.inter_size, 1.2, 30u);
    const fourier plant_fitness(sp.p_gene_size, sp.inter_size, 1.2, 30u);

    const unsigned ancestor_size = 50;
    auto [b_pop, p_pop] = random_haplotypes(sp, bacteria_fitness, plant_fitness, ancestor_size, 1.25, 1.4);
    b_pop.reserve(sp.b_pop_size);
    p_pop.reserve(sp.p_pop_size);

    std::experimental::filesystem::create_directory(sp.out_path);

    std::mt19937 random_generator((std::random_device())());

    std::cout << sp.p_pop_size << std::endl;

    population bacteria(sp.bm, sp.bv, bacteria_fitness, sp.b_pop_size, ancestor_size, 0, b_pop);
    population plant(sp.bm, sp.bv, plant_fitness, sp.p_pop_size, ancestor_size, 0, p_pop);

    simulation sim{ bacteria, plant, sp.generations, sp.out_path };

    sim.run();
}

int main(int argc, char *argv[])
{
    std::ios_base::sync_with_stdio(false);

    std::vector<std::tuple<unsigned, unsigned>> pop_sizes({ { 100000, 1000 }, { 100000, 100 } });
    std::vector<std::tuple<double, double>> p_mut_params({ { 0.06, 0.15 } });
    std::vector<std::tuple<double, double>> b_mut_params({ { 0.005, 0.05 }, { 0.005, 0.015 }, { 0.005, 0.12 } });
    std::vector<std::tuple<unsigned, unsigned>> fitness_params({ { 5, 2 } });

    std::vector<simulation_params> sp;

    if (argc != 2) {
        std::cout << "try again you corky dork" << std::endl;
        return 0;
    }

    std::experimental::filesystem::path out_path(std::experimental::filesystem::current_path());
    out_path.append("mut_prob_dif_pic");
    std::experimental::filesystem::create_directory(out_path);

    for (auto id = 0; id < 4000; id++)
        for (auto& ps : pop_sizes)
            for (auto& b_params : b_mut_params)
                for (auto& p_params : p_mut_params)
                    for (auto& f_params : fitness_params) {
                        std::experimental::filesystem::path out_path_in(out_path);

                        std::ostringstream oss;
                        oss << "(" << std::get<0>(f_params) << "," << std::get<0>(f_params) << "," << std::get<1>(f_params) << "),";
                        print(oss, ps);
                        oss << ",";
                        print(oss, b_params);
                        oss << ",";
                        print(oss, p_params);
                        oss << ",2500," << (std::random_device())();
                        out_path_in.append(oss.str());

                        sp.push_back({ 2500u, std::get<0>(f_params), std::get<0>(f_params), std::get<1>(f_params),
                                       std::get<0>(ps), std::get<1>(ps),
                                       std::get<0>(b_params) / std::get<0>(f_params), std::get<1>(b_params),
                                       std::get<0>(p_params) / std::get<0>(f_params), std::get<1>(p_params), out_path_in });
                }


    ctpl::thread_pool pool(std::strtol(argv[1], nullptr, 10));

    for (auto& i : sp)
        pool.push(run_simulation, i);

    return 0;
}