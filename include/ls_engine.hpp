#pragma once
#include "global.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include "stats_stream.hpp"
#include <vector>
#include <unordered_set>
#include <deque>

class ls_engine
{
    using string = std::string;

public:
    // ======= Meta data  ============
    struct move
    {
        int op1, op2;
        bool fw;

        string to_string(instance *inst);

        bool check(solution &s);
    };

    class move_mem_data
    {
        vi block_idxs;
        static std::hash<int> hsh_int;

    public:
        bool fw;
        size_t seed = 0;

        void init(bool _fw);

        struct eqp_cmp
        {
            bool operator()(const move_mem_data &move_mem1, const move_mem_data &move_mem2) const;
        };

        void add_value(int x);

        struct move_hash
        {
            size_t operator()(const move_mem_data &move) const;
        };

        string to_string();
    };

    // ====== Memory data =========
    struct ls_stats
    {
        stats_stream<long> one_round_time_stats = stats_stream<long>("round time [ms]");
        long total_time = 0, n_ls_moves = 0;
        std::vector<move> moves_memory;
        long n_moves_no_local_imp = 0;
        long n_moves_no_global_imp = 0;

        vvi best_sol_found; //machine perm representation
        int best_fo_found = I32_INF;
        bool stop_because_no_move = false;

        std::unordered_set<move_mem_data, move_mem_data::move_hash, move_mem_data::eqp_cmp> tabu_set;
        std::deque<move_mem_data> tabu_list;

        void reset(solution &s);

    } stats;

    // ===== LS parameters ================
    struct ls_parameters
    {
        bool debug_mode = false;

        //Neighborhood Exploration
        bool choose_rand_first_block = true;
        bool explore_all_blocks = true;
        bool stop_first_improve = false;
        bool not_down_moves = false;
        bool break_even_moves = false;

        //Stop Criteria
        bool restrict_n_iterations = false;
        int max_n_moves = 10000;
        bool restrict_max_no_global_improve = true;
        int max_moves_no_global_improve = 1000;

        //Tabu Search
        int tabu_list_size = 15;
        bool change_tabu_list_size = false;
        int tabu_list_percent_change = 20; //[%]
        bool strong_tabu_rule = false;
        double tabu_list_factor = 2;

        void compute(solution &s);

        ls_parameters();

    } params;

private:
public:
    
    ls_engine();

    //Parameters Input
    void config(ls_parameters _params = ls_parameters());

    //String signature
    string signature();

    void run_ls(solution &initial_sol);

    int best_fo_found();

    //returns the job permutation in each machine
    vvi best_sol_found();

private:
    // ======== Locas Search Algorithms =============
    int n, m, nop;
    solution *s;
    vi hb_dist, tb_dist; //block distances
    std::vector<bool> mark;   //aux mark vector
    std::vector<int> aux_buf;
    vi explored_blocks;

    bool tabu_acc(move_mem_data &move_tabu_data, int eval);

    void tabu_add(move_mem_data &move_tabu_data);

    void create_tabu_move_data(move &move, vi &ms, vi &mp,
                               std::vector<move_mem_data> &v_tabu_move_data);

    void tabu_list_size_update();

    bool test_acyclic_move(int u, int v, bool is_fw);

    void explore_block_moves(vi &block_idxs,
                             vi &js, vi &jp, vi &ms, vi &mp, vi &head_dist, vi &tail_dist, vi &time,
                             int &best_eval, move &best_move,
                             bool fw);

    bool search_block(int block_begin, move &best_move, int &best_eval);

    bool search_move(move &selected_move);

    void update_tsv(int u, int v, bool fw,
                    int &u_init_pos, int &v_init_pos);

    void apply_move(move &move_apply);

    bool stop_criterion_test();

    void keep_best_sol_test();

    void ls_start(solution &s);
};