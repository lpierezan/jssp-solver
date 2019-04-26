#pragma once
#include "global.hpp"
#include "instance.hpp"
#include "ls_engine.hpp"
#include "solution.hpp"
#include "solution_generator.hpp"
#include <chrono>
#include <vector>

class mh_engine
{
    instance *inst;

public:
    // ==== MH Parameters ======
    struct mh_params
    {
        ls_engine ls;
        bool log_mode = false;

        //==== Stop Criteria
        bool stop_max_no_imp_it = true;
        int max_no_imp_it_limit = 100;
        bool stop_max_it = false;
        int max_it_limit = 1000;
        bool stop_time_limit = false;
        long time_limit_sec = 60;
        //===========

        // ======== Genetic Population ===========
        bool improve_after_generation = true;
        int max_parents_no_offspring = 15;
        int initial_pop_size = 100;
        int curr_pop_size = 100;
        int n_off_to_gen = 70;
        int n_elite = 20;
        bool second_parent_non_elite = true;
        ls_engine::ls_parameters ls_params_indv_gen;
        solution_generator sol_generator;
        // ===============================

        // ======== Path relink ========
        bool generate_both_directions = true;
        double break_path_blocks = 15;
        bool cut_pathrelink_borders = true;
        double denominator_cut_pathrelink_borders = 6;
        bool repair_infeasible = true;

        bool do_offspring_intensive_ls = true;
        bool do_path_eval_ls = true;
        ls_engine::ls_parameters ls_params_path_eval;
        ls_engine::ls_parameters ls_params_offspring;
        //=============================

        mh_params();

        void init(instance &inst);
    } params;

    // ==== MH Memory/Stats Data
    struct mh_stats
    {
        //Solutions data
        int best_fo_found = I32_INF;
        solution best_sol_found; //only set at the end

        //Time measures
        long total_time = 0;
        std::chrono::time_point<std::chrono::steady_clock> start_time;
        stats_stream<long> it_time_stream;

        int curr_iteration = 0, n_iter_no_global_improve = 0;
        stats_stream<long long> parents_dist;
        int n_infeasible = 0;

        void init(instance &inst);

    } stats;

    // === Meta data ========
    class individual
    {
    public:
        vvi mch_perm;
        int fitness;

        individual(solution &s);

        std::string to_string(instance &inst);

        static individual make_invalid();

    private:
        individual();
    };

    class population
    {
        std::vector<individual *> indv_pointers;
        int best_fit = I32_INF;
        bool sorted = false;

    public:
        ~population();

        int best_fitness();

        int size();

        void clear();

        void add(const individual &_indv);

        //empty offspring
        void consume(population &off);

        void pop_back();

        void sort();

        void keep_best(int n_best);

        const individual *select_bests(int n_elite);

        const individual *select_worst(int n_worse);

        std::string to_string(instance &inst);

        individual get_best();

    private:
        void clear_no_delete();
    };

public:
    void config(mh_params _mh_params = mh_params());

    solution run(instance &inst_);

    long get_best_fo();

    solution get_best_sol();

    std::string signature();

private:
    vvi aux1; //[m,n] pre-allocated

    individual gen_new_individual();

    void gen_init_pop(population &pop);

    void update_best_fitness(population &pop);

    bool stop_mh(population &pop);

    //must return <elite,normal>
    std::pair<const individual *, const individual *> select_parents(population &pop);

    void apply_path_move(int &curr_dist, vvi &curr_indv, const vvi &target_indv, std::vector<std::pair<int, int>> &wrong_pos, vvi &f);

    int distance(const individual &p1, const individual &p2);

    void generate_childrens_path_relink(const individual &p1, const individual &p2, population &off);

    void generate_childrens(const individual &elite, const individual &normal, population &off);

    void generate_offspring(population &pop, population &off);

    void update_pop(population &pop, population &off);

    void run_mh();
};