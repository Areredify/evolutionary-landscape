#pragma once

#include <random>
#include <set>
#include <vector>
#include <deque>
#include <sstream>
#include <unordered_map>
#include <experimental/filesystem>

double rnd0(std::mt19937 &random_generator);
double rnd1(std::mt19937 &random_generator);
double nrm0(std::mt19937 &random_generator);
double sqr(double x);

double distance(const std::vector<double> & v1, 
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

    fourier(int input_size, int output_size, double power, int func_num);

    std::vector<double> operator()(const std::vector<double> & gene) const override;

    fourier(const fourier & f) = default;
    fourier(fourier && f) = default;
    fourier& operator=(const fourier& f) = delete;
    fourier& operator=(fourier&& f) = delete;

    virtual ~fourier() = default;
};

struct entry {
    std::vector<double> gene, t_gene;
};

typedef std::unordered_map<size_t, entry> database;

class organism {
public:
    unsigned id;
    std::deque<unsigned> ancestors;

    organism(const unsigned id, std::deque<unsigned> ancestors) : id(id), ancestors(std::move(ancestors)) {
    }
};

class population {
public:
    const double m, v;
    const fitness& fitness_func;
    const unsigned size;
    const unsigned ancestor_size;
    database db;
    unsigned db_last_id;
    std::vector<organism> pop;

    void cleanup_database() {
        database new_db;

        for (auto & i : pop) {
            if (new_db.count(i.id) == 0)
                new_db[i.id] = db[i.id];
            for (auto & j : i.ancestors)
                if (new_db.count(j) == 0)
                    new_db[j] = db[j];
        }

        db = std::move(new_db);
    }

    auto min_max() {
        auto mn = db[pop[0].id].t_gene, mx = db[pop[0].id].t_gene;

        for (auto & i : pop) {
            auto& t_gene = db[i.id].t_gene;
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
            auto& t_gene = db[i.id].t_gene;
            for (auto j = 0u; j < fitness_func.output_size; j++) {
                current_index = (current_index * power + static_cast<unsigned>((t_gene[j] - mn[j] - 1e-6) / (mx[j] - mn[j]) * power));
            }
            discrete_representation[current_index] += 1.0 / pop.size();
        }

        return discrete_representation;
    }

    std::string stats() {
        double speed = 0;

        for (auto & i : pop)
            if (!i.ancestors.empty())
                speed += distance(db[i.id].gene, db[i.ancestors.front()].gene) / i.ancestors.size();

        std::ostringstream oss;

        oss << speed / pop.size() / v;

        return oss.str();
    }

    friend std::ostream& operator<<(std::ostream & stream, const population& pop);

    organism make_offspring(const organism& o, std::mt19937 &random_generator) {
        auto gene = db[o.id].gene;

        auto flag = false;
        for (auto &i : gene)
            if (rnd0(random_generator) < m) {
                flag = true;
                i += nrm0(random_generator) * v;
            }

        auto res = o;
        res.ancestors.push_back(o.id);
        if (res.ancestors.size() > ancestor_size)
            res.ancestors.pop_front();

        if (!flag)
            return res;

        db[db_last_id] = entry{ gene, fitness_func(gene) };
        res.id = db_last_id;
        db_last_id++;

        return res;
    }
};

class simulation {
public:
    population &bacteria_pop, &plant_pop;
    unsigned generations;

    std::experimental::filesystem::path out_path;

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
