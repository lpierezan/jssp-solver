#pragma once
#include "global.hpp"
#include <string>
#include <random>
#include <vector>
#include "operation.hpp"

class instance
{
    using string = std::string;

    string inst_filepath;
    std::default_random_engine my_rand_gen = std::default_random_engine(0);

public:
    int n, m, nop;
    vvi time, machine; //[job idx][pos idx]
    vvi jm_to_pos;

    vi time_idx, mch_idx, js_idx, jp_idx, job_idx;

    void initialize(string filepath);

    // ================== Algorithms =============================

    void build_job_neighborhood();

    void build_mch_neighborhood(vvi &mch_perm, vi &ms, vi &mp);

    bool dfs(operation op, std::vector<std::vector<char>> &mark, vi &ms, std::vector<operation> &tsv);
    
    bool compute_fo(vvi mch_perm, int &fo);

    //Build active scheduling given operation order. O(n*m*log(n))
    void build_active_scheduling(vi &op_perm, vvi &mch_perm);

    // ================== Utilities =============================
    void print();

    string mch_perm_to_str(const vvi &mch_perm);

    inline std::pair<int, int> pair_job_pos(int idx)
    {
        int pi = idx % m;
        if (pi == 0)
            pi = m;
        int ji = ((idx - pi) / m) + 1;
        return std::make_pair(ji, pi);
    }

    inline std::pair<int, int> pair_job_mch(int idx)
    {
        return std::make_pair(job_idx[idx], mch_idx[idx]);
    }

    inline int idx(int ji, int pi)
    {
        return (ji - 1) * m + pi;
    }

    inline int idx(operation op)
    {
        return idx(op.ji, op.pi);
    }

    bool hard_check(vvi const &machine_perm, int cmp_fo = I32_INF);

    void job_perm_to_op_perm(vi &job_perm);

    void restart_pi_gen(int seed = 0);

    vi generate_random_pi(int min, int max);

    string op_vector_to_str(vi &op_vec, string sep = " ");

    bool check_mch_ng(int v, vi &ms, vi &mp);
};