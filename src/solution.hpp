#pragma once
#include "global.hpp"
#include "instance.hpp"
#include <vector>
#include <utility>

class solution
{
    
public:
    instance *my_inst = 0;
    int n, m, nop;
    vi mp_idx, ms_idx, head_dist, tail_dist; //[idx] => ()
    vi tsv, tsv_inv;
    vvi my_mch_perm;
    int fo_online = I32_INF;
    bool valid = false;

    solution();
    solution(instance &inst);

    void associate(instance &inst);

    bool config(vvi &mch_perm);

    bool config(vi &op_perm);

    bool repair();

    int get_fo();

    // ============ Utilities ======================
    int inline mch_idx(int op) { return my_inst->mch_idx[op]; }

    std::pair<int, int> get_succ(int op);

    std::pair<int, int> get_pred(int op);

    int inline js(int idx) { return my_inst->js_idx[idx]; }

    int inline jp(int idx) { return my_inst->jp_idx[idx]; }

    int inline idx(int ji, int pi) { return my_inst->idx(ji, pi);}

    bool inline is_critical(int op)
    {
        //0 -> op + op -> end
        return ((tail_dist[op] + head_dist[op] - my_inst->time_idx[op]) == fo_online);
    }

    bool inline is_block_start(int op)
    {
        return (is_critical(op) && (mp_idx[op] == 0 || !is_critical(mp_idx[op])));
    }

    void print(bool print_mach_perm = false);

private:
    //Dfs used for topological sort
    bool ts_dfs(int op_idx);

    //Build topological sort vector
    std::vector<char> mark;
    bool build_ts();

    //Build distance vectors from head and tail
    bool build_dists();
};