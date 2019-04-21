#pragma once

#include <random>
#include <set>
#include <vector>
#include <deque>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <array>

#include "utils.h"

double rnd0(std::mt19937 &random_generator);
double rnd1(std::mt19937 &random_generator);
double nrm0(std::mt19937 &random_generator);
double sqr(const double x);

const auto max_gene_size = 10;

using gene_c = truncated_array<double, max_gene_size>;

double distance(const gene_c &v1,
                const gene_c &v2);
double distance2(const gene_c &v1,
                 const gene_c &v2);
gene_c rnd1_vec(unsigned len, std::mt19937 &random_generator);

template <class TupType, size_t... I>
void print(std::ostream &os, TupType &_tup, std::index_sequence<I...>)
{
    os << "(";
    (..., (os << (I == 0 ? "" : ",") << std::get<I>(_tup)));
    os << ")";
}

template <class... T>
void print(std::ostream &os, const std::tuple<T...> &_tup)
{
    print(os, _tup, std::make_index_sequence<sizeof...(T)>());
}

class fitness
{
  public:
    const unsigned input_size, output_size;
    virtual gene_c operator()(const gene_c &gene) const = 0;

    fitness(const unsigned input_size, const unsigned output_size)
        : input_size(input_size), output_size(output_size)
    {
    }

    fitness(const fitness &f) = default;
    fitness(fitness &&f) = default;
    fitness &operator=(const fitness &f) = delete;
    fitness &operator=(fitness &&f) = delete;

    virtual ~fitness() = default;
};

class fourier : public fitness
{
  public:
    struct cell
    {
        double a, k, b;
    };

    const unsigned func_num;
    std::vector<cell> transformer;

    unsigned index(int o, int i, int f) const;

    fourier(const unsigned input_size, const unsigned output_size, const double power, const unsigned func_num);

    gene_c operator()(const gene_c &gene) const override;

    fourier(const fourier &f) = default;
    fourier(fourier &&f) = default;
    fourier &operator=(const fourier &f) = delete;
    fourier &operator=(fourier &&f) = delete;

    virtual ~fourier() = default;
};

const auto max_ancestor_size = 10;
class organism
{
  public:
    gene_c gene;
    gene_c t_gene;
    truncated_array<organism *, max_ancestor_size> ancestors;

    organism *k_ancestor;

    organism() = default;
    //explicit organism(gene_c gene) : gene(gene) {
    //}
};

//typedef std::vector<std::shared_ptr<organism>> org_vec;

class population
{
  public:
    const double m, v;
    const fitness &fitness_func;
    const unsigned size;
    const unsigned ancestor_size;
    unsigned speed_scale;
    std::vector<organism> pop_storage;
    unsigned cur_pop_row;
    organism *cur_pop;
    std::vector<double> speed_reference;
    unsigned *epoch;

    population(const double m, const double v, const fitness &fitness_func,
               const unsigned size, const unsigned ancestor_size,
               const unsigned speed_scale, std::vector<organism> pop_storage) : m(m), v(v), fitness_func(fitness_func), size(size), ancestor_size(ancestor_size), speed_scale(speed_scale),
                                                                                pop_storage(std::move(pop_storage)), cur_pop_row(0), speed_reference(ancestor_size), epoch(nullptr)
    {
        this->pop_storage.resize(size * (ancestor_size + 1), {{this->pop_storage[0].gene.size}, {this->pop_storage[0].t_gene.size}, this->pop_storage[0].ancestors.size, nullptr});
        cur_pop = get_cur_pop();

        std::mt19937 gen(std::random_device{}());
        const auto iters = 10000u;

        for (auto i = 0u; i < iters; i++)
        {
            double sum = 0;
            for (auto j = 0u; j < ancestor_size; j++)
            {
                if (rnd0(gen) < m)
                    sum += nrm0(gen) * v;
                speed_reference[j] += sqr(sum);
            }
        }

        for (auto &i : speed_reference)
            i = i / iters * fitness_func.input_size;
    }

    organism *get_cur_pop()
    {
        return pop_storage.data() + size * cur_pop_row;
    }

    void advance_pop()
    {
        cur_pop_row = (cur_pop_row + 1) % (ancestor_size + 1);
        cur_pop = get_cur_pop();
    }

    auto min_max()
    {
        const auto cur_pop = get_cur_pop();

        auto mn = cur_pop[0].t_gene, mx = cur_pop[0].t_gene;

        for (auto i = 1u; i < size; i++)
        {
            auto &t_gene = cur_pop[i].t_gene;
            for (auto j = 0u; j < fitness_func.output_size; j++)
            {
                mn[j] = std::min(mn[j], t_gene[j]);
                mx[j] = std::max(mx[j], t_gene[j]);
            }
        }

        return std::tuple{mn, mx};
    }

    std::vector<double> labeling(const gene_c &mn, const gene_c &mx) const
    {
        const unsigned power = 15;
        std::vector<double> discrete_representation(static_cast<unsigned>(pow(power, fitness_func.output_size) + 0.1));

        for (auto i = 0u; i < size; i++)
        {
            unsigned current_index = 0;
            auto &t_gene = cur_pop[i].t_gene;
            for (auto j = 0u; j < fitness_func.output_size; j++)
            {
                auto x = static_cast<unsigned>((t_gene[j] - mn[j] - 1e-6) / (mx[j] - mn[j]) * power);
                current_index = current_index * power + x;
            }
            discrete_representation[current_index] += 1.0 / size;
        }

        return discrete_representation;
    }

    double pop_speed()
    {
        double speed = 0;

        for (auto i = 0u; i < size; i++)
            if (cur_pop[i].k_ancestor != nullptr)
                speed += distance2(cur_pop[i].gene, cur_pop[i].k_ancestor->gene);

        return speed / size / speed_reference[std::min(ancestor_size, speed_scale) - 1];
    }

    double pop_dispersion()
    {
        gene_c center{fitness_func.input_size};

        for (auto i = 0u; i < size; i++)
        {
            for (auto j = 0u; j < fitness_func.input_size; j++)
                center[j] += cur_pop[i].gene[j];
        }

        for (auto i = 0u; i < fitness_func.input_size; i++)
            center[i] /= size;

        double res = 0;

        for (auto i = 0u; i < size; i++)
            res += distance2(center, cur_pop[i].gene);

        return res / size;
    }

    std::string stats()
    {
        std::ostringstream oss;

        oss << pop_speed() << " " << pop_dispersion();

        return oss.str();
    }

    friend std::ostream &operator<<(std::ostream &stream, const population &pop);

    void make_offspring(organism *o, organism *output, std::mt19937 &random_generator) const
    {

        auto flag = false;
        for (auto i = 0u; i < output->gene.size; i++)
            if (rnd0(random_generator) < m)
            {
                flag = true;
                output->gene[i] += nrm0(random_generator) * v;
            }

        if (flag)
            output->t_gene = fitness_func(output->gene);
        else
            output->t_gene = o->t_gene;

        unsigned cur_ind = 0;
        auto cur_ptr = o;

        for (auto i = 0u; i < output->ancestors.size; i++)
        {
            if (cur_ptr == nullptr)
                break;

            output->ancestors[i] = cur_ptr;
            cur_ptr = cur_ptr->ancestors[cur_ind++];
        }

        cur_ind = 0;
        cur_ptr = output;

        for (auto k = std::min(ancestor_size, *epoch); k > 0; k /= 2, cur_ind++)
            if (k & 1)
                cur_ptr = cur_ptr->ancestors[cur_ind];
        output->k_ancestor = cur_ptr;
    }
};

class simulation
{
  public:
    population &bacteria_pop, &plant_pop;
    unsigned generations;

    std::filesystem::path out_path;

    double get_similarity() const;
    double get_distance(std::mt19937 &random_generator) const;

    static double get_prob(double mean);
    void evolve(std::mt19937 &random_generator);
    void run();
};

template <class InputIt>
std::vector<InputIt> KMax(unsigned k, InputIt first, InputIt last)
{
    class comparator
    {
      public:
        bool operator()(const InputIt &first, const InputIt &last) const
        {
            return *first < *last;
        }
    } foo;

    std::multiset<InputIt, comparator> res(foo);
    for (InputIt i = first; i != last; ++i)
    {
        res.insert(i);
        if (res.size() > k)
            res.erase(res.begin());
    }
    return std::vector<InputIt>(res.begin(), res.end());
}

template <class InputIt>
std::vector<InputIt> KMin(unsigned k, InputIt first, InputIt last)
{
    class comparator
    {
      public:
        bool operator()(const InputIt &first, const InputIt &last) const
        {
            return *first > *last;
        }
    } foo;

    std::multiset<InputIt, comparator> res(foo);
    for (auto i = first; i != last; ++i)
    {
        res.insert(i);
        if (res.size() > k)
            res.erase(res.begin());
    }

    return std::vector<InputIt>(res.begin(), res.end());
}

class simulation_params
{
  public:
    unsigned generations;
    unsigned b_gene_size, p_gene_size;
    unsigned inter_size;
    unsigned b_pop_size, p_pop_size;

    double bv, bm, pv, pm;

    std::filesystem::path out_path;
};

template <class Func>
auto itemwise_func(const gene_c &v1, const gene_c &v2, Func f)
{
    gene_c res;
    res.size = v1.size;

    for (auto i = 0u; i < v1.size; i++)
        res[i] = f(v1[i], v2[i]);

    return res;
}
