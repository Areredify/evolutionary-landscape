#pragma once

#include <random>
#include <set>
#include <vector>
#include <deque>
#include <sstream>
#include <unordered_map>
#include <experimental/filesystem>
#include <cassert>

double rnd0(std::mt19937 &random_generator);
double rnd1(std::mt19937 &random_generator);
double nrm0(std::mt19937 &random_generator);
std::vector<double> rnd1_vec(unsigned len, std::mt19937 &random_generator);
double sqr(const double x);

double distance(const std::vector<double> & v1,
                const std::vector<double> & v2);
double distance2(const std::vector<double> & v1,
                 const std::vector<double> & v2);

void print_gene(const std::vector<double> & v);

template<class TupType, size_t... I>
void print(std::ostream & os, TupType& _tup, std::index_sequence<I...>)
{
    os << "(";
    (..., (os << (I == 0 ? "" : ",") << std::get<I>(_tup)));
    os << ")";
}

template<class... T>
void print(std::ostream & os, const std::tuple<T...>& _tup)
{
    print(os, _tup, std::make_index_sequence<sizeof...(T)>());
}

class fitness {
public:
    const unsigned input_size, output_size;
    virtual std::vector<double> operator()(const std::vector<double> & gene) const = 0;

    fitness(const unsigned input_size, const unsigned output_size)
        : input_size(input_size), output_size(output_size) {
    }

    fitness(const fitness & f) = default;
    fitness(fitness && f) = default;
    fitness& operator=(const fitness& f) = delete;
    fitness& operator=(fitness&& f) = delete;

    virtual ~fitness() = default;
};

class fourier : public fitness {
public:
    struct cell {
        double a, k, b;
    };

    const unsigned func_num;
    std::vector<cell> transformer;

    unsigned index(int o, int i, int f) const;

    fourier(const unsigned input_size, const unsigned output_size, const double power, const unsigned func_num);

    std::vector<double> operator()(const std::vector<double> & gene) const override;

    fourier(const fourier & f) = default;
    fourier(fourier && f) = default;
    fourier& operator=(const fourier& f) = delete;
    fourier& operator=(fourier&& f) = delete;

    virtual ~fourier() = default;
};

class organism {
public:
    std::vector<double> gene;
    std::vector<double> t_gene;
    std::vector<std::shared_ptr<organism>> ancestors;

    std::shared_ptr<organism> k_ancestor;

    organism() = default;
    explicit organism(std::vector<double> gene) : gene(std::move(gene)) {
    }
};

typedef std::vector<std::shared_ptr<organism>> org_vec;

class population {
public:
    const double m, v;
    const fitness& fitness_func;
    const unsigned size;
    const unsigned ancestor_size;
    unsigned speed_scale;
    std::vector<std::shared_ptr<organism>> pop;
    std::vector<double> speed_reference;
    unsigned * epoch;

    population(const double m, const double v, const fitness & fitness_func,
               const unsigned size, const unsigned ancestor_size, 
               const unsigned speed_scale, std::vector<std::shared_ptr<organism>> pop) :
    m(m), v(v), fitness_func(fitness_func), size(size), ancestor_size(ancestor_size), speed_scale(speed_scale),
    pop(std::move(pop)), speed_reference(ancestor_size), epoch(nullptr) {
        std::mt19937 gen(std::random_device{}());
        const auto iters = 1000000u;

        for (auto i = 0u; i < iters; i++) {
            double sum = 0;
            for (auto j = 0u; j < ancestor_size; j++) {
                if (rnd0(gen) < m)
                    sum += nrm0(gen) * v;
                speed_reference[j] += sqr(sum);
            }
        }

        for (auto &i : speed_reference)
            i = i / iters * fitness_func.input_size;
    }

    auto min_max() {
        auto mn = pop[0]->t_gene, mx = pop[0]->t_gene;

        for (auto & i : pop) {
            auto& t_gene = i->t_gene;
            for (auto j = 0u; j < fitness_func.output_size; j++) {
                mn[j] = std::min(mn[j], t_gene[j]);
                mx[j] = std::max(mx[j], t_gene[j]);
            }
        }

        return std::tuple{ mn, mx };
    }

    std::vector<double> labeling( const std::vector<double>& mn, const std::vector<double>& mx) {
        const unsigned power = 15;
        std::vector<double> discrete_representation(static_cast<unsigned>(pow(power, fitness_func.output_size) + 0.1));

        for (auto & i : pop) {
            unsigned current_index = 0;
            auto& t_gene = i->t_gene;
            for (auto j = 0u; j < fitness_func.output_size; j++) {
                auto x = static_cast<unsigned>((t_gene[j] - mn[j] - 1e-6) / (mx[j] - mn[j]) * power);
                current_index = current_index * power + x;
            }
            discrete_representation[current_index] += 1.0 / pop.size();
        }

        return discrete_representation;
    }

    double pop_speed() {
        double speed = 0;

        for (auto & i : pop)
            if (i->k_ancestor != nullptr)
                speed += distance2(i->gene, i->k_ancestor->gene);

        return speed / pop.size() / speed_reference[std::min(ancestor_size, speed_scale) - 1];
    }

    double pop_dispersion() {
        std::vector<double> center(fitness_func.input_size);

        for (auto & i : pop) {
            for (auto j = 0u; j < fitness_func.input_size; j++)
                center[j] += i->gene[j];
        }

        for (auto & i : center)
            i /= size;

        double res = 0;

        for (auto & i : pop)
            res += distance2(center, i->gene);

        return res / size;
    }

    std::string stats() {
        

        std::ostringstream oss;

        oss << pop_speed() << " " << pop_dispersion();

        return oss.str();
    }

    friend std::ostream& operator<<(std::ostream & stream, const population& pop);

    std::shared_ptr<organism> make_offspring(const std::shared_ptr<organism> o, std::mt19937 &random_generator) {
        auto res = std::make_shared<organism>(o->gene);

        auto flag = false;
        for (auto &i : res->gene)
            if (rnd0(random_generator) < m) {
                flag = true;
                i += nrm0(random_generator) * v;
            }

        if (flag)
            res->t_gene = fitness_func(res->gene);
        else
            res->t_gene = o->t_gene;

        unsigned cur_ind = 0;
        auto cur_ptr = o;
        res->ancestors.resize(o->ancestors.size());

        for (auto & i : res->ancestors) {
            if (cur_ptr == nullptr)
                break;

            i = cur_ptr;
            cur_ptr = cur_ptr->ancestors[cur_ind++];
        }

        cur_ind = 0;
        cur_ptr = res;

        for (auto k = std::min(ancestor_size, *epoch); k > 0; k /= 2, cur_ind++)
            if (k & 1)
                cur_ptr = cur_ptr->ancestors[cur_ind];
        res->k_ancestor = cur_ptr;

        return res;
    }
};

class simulation {
public:
    population &bacteria_pop, &plant_pop;
    unsigned generations;

    std::experimental::filesystem::path out_path;

    double get_similarity() const;
    double get_distance(std::mt19937 & random_generator) const;

    static double get_prob(double mean);
    void evolve(std::mt19937 &random_generator);
    void run();
};

template< class InputIt >
std::vector<InputIt> KMax(unsigned k, InputIt first, InputIt last)
{
    class comparator {
    public:
      bool operator()(const InputIt & first, const InputIt & last) const {
          return *first < *last;
      }
    } foo;
    
    std::multiset<InputIt, comparator> res(foo);
    for (InputIt i = first; i != last; ++i) {
        res.insert(i);
        if (res.size() > k)
            res.erase(res.begin());
    } 
    return std::vector<InputIt>(res.begin(), res.end());
}

template< class InputIt >
std::vector<InputIt> KMin(unsigned k, InputIt first, InputIt last)
{
    class comparator {
    public:
        bool operator()(const InputIt & first, const InputIt & last) const {
            return *first > *last;
        }
    } foo;

    std::multiset<InputIt, comparator> res(foo);
    for (auto i = first; i != last; ++i) {
      res.insert(i);
      if (res.size() > k)
        res.erase(res.begin());
    }

    return std::vector<InputIt>(res.begin(), res.end());
}

class simulation_params {
public:
    unsigned generations;
    unsigned b_gene_size, p_gene_size;
    unsigned inter_size;
    unsigned b_pop_size, p_pop_size;

    double bv, bm, pv, pm;

    std::experimental::filesystem::path out_path;
};

template <class T, class Func>
auto itemwise_func(const std::vector<T> & v1, const std::vector<T> & v2, Func f) {
    std::vector<T> res(v1.size());

    for (auto i = 0u; i < v1.size(); i++)
        res[i] = f(v1[i], v2[i]);

    return res;
}
