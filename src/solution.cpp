#include "global.hpp"
#include "solution.hpp"
#include <iostream>
#include <utility>
#include <algorithm>

using namespace std;

solution::solution()
{
    valid = false;
    my_inst = 0;
};

solution::solution(instance &inst)
{
    associate(inst);
}

void solution::associate(instance &inst)
{
    my_inst = &inst;
    n = my_inst->n;
    m = my_inst->m;
    nop = my_inst->nop;
    valid = false;
}

bool solution::config(vvi &mch_perm)
{
    if (my_inst == 0)
    {
        cout << "No instance associated error.";
        exit(1);
    }
    fo_online = I32_INF;
    my_mch_perm = mch_perm;

    //building machine pred. and succ.
    my_inst->build_mch_neighborhood(mch_perm, ms_idx, mp_idx);

    //build distances
    auto ok = build_dists();
    valid = ok;

    return ok;
}

bool solution::config(vi &op_perm)
{
    my_mch_perm.clear();
    my_inst->build_active_scheduling(op_perm, my_mch_perm);
    return config(my_mch_perm);
}

bool solution::repair()
{
    if (valid)
        return true; //nothing todo

    vi job_perm;
    job_perm.reserve(nop);
    int mi, pi;
    for (pi = 0; pi < n; pi++)
    {
        for (mi = 0; mi < m; mi++)
        {
            int ji = my_mch_perm[mi][pi];
            job_perm.push_back(ji);
        }
    }
    my_inst->job_perm_to_op_perm(job_perm);
    bool ok = config(job_perm);
    if (!ok)
    {
        cout << "not expected";
        exit(1);
    }
    return ok;
}

int solution::get_fo()
{
    return fo_online;
}

pair<int, int> solution::get_succ(int op)
{
    return make_pair(my_inst->js_idx[op], ms_idx[op]);
}

pair<int, int> solution::get_pred(int op)
{
    return make_pair(my_inst->jp_idx[op], mp_idx[op]);
}

void solution::print(bool print_mach_perm)
{
    if (!valid)
    {
        cout << "Invalid solution" << endl;
        return;
    }

    cout << "Makespan = " << fo_online << endl;
    if (print_mach_perm)
    {
        cout << my_inst->mch_perm_to_str(my_mch_perm);
    }
}

//Dfs used for topological sort
bool solution::ts_dfs(int op_idx)
{
    mark[op_idx] = 1;
    bool ok;
    int succ;

    auto succ_pair = get_succ(op_idx);
    for (int ns = 1; ns <= 2; ns++)
    {
        if (ns == 1)
        {
            //Job succesor
            succ = succ_pair.first;
        }
        else
        {
            //Machine sucessor
            succ = succ_pair.second;
        }

        //Process sucessor
        if (succ != 0)
        {
            if (!mark[succ])
            {
                ok = ts_dfs(succ);
                if (!ok)
                    return false;
            }
            else if (mark[succ] == 1)
            {
                return false;
            }
            else
            {
                //mark[succ] == 2 nothing to do
            }
        }
    }

    mark[op_idx] = 2;
    tsv.push_back(op_idx);
    return true;
}

//Build topological sort vector

bool solution::build_ts()
{
    bool ok;

    tsv.clear();
    if (mark.size() == 0)
    {
        mark.resize(nop + 1, 0);
    }
    else
    {
        fill(ALL(mark), 0);
    }

    int op_idx;
    for (int ji = 1; ji <= n; ji++)
    {
        op_idx = my_inst->idx(ji, 1);
        if (!mark[op_idx])
        {
            ok = ts_dfs(op_idx);
            if (!ok)
                return false;
        }
    }

    if (tsv.size() != nop)
        return false;

    reverse(ALL(tsv));
    tsv_inv.resize(tsv.size() + 1);
    for (int i = 0; i < tsv.size(); i++)
        tsv_inv[tsv[i]] = i;

    //ok
    return true;
}

//Build distance vectors from head and tail
bool solution::build_dists()
{
    auto ok = build_ts();
    if (!ok)
        return false;

    tail_dist.resize(nop + 1, 0);
    head_dist.resize(nop + 1, 0);
    head_dist[0] = 0;

    //Head distances
    for (int i = 0; i < tsv.size(); i++)
    {
        int op = tsv[i], pred;
        head_dist[op] = 0;

        auto pred_pair = get_pred(op);
        for (int nv = 1; nv <= 2; nv++)
        {
            if (nv == 1)
                pred = pred_pair.first;
            else
                pred = pred_pair.second;

            if (pred != 0)
                head_dist[op] = max(head_dist[op], head_dist[pred]);
        }
        head_dist[op] += my_inst->time_idx[op];
    }

    //Tail distances
    tail_dist[0] = 0;
    fo_online = 0;
    for (int i = tsv.size() - 1; i >= 0; i--)
    {
        int op = tsv[i], succ;
        tail_dist[op] = 0;

        auto succ_pair = get_succ(op);
        for (int nv = 1; nv <= 2; nv++)
        {
            if (nv == 1)
                succ = succ_pair.first;
            else
                succ = succ_pair.second;

            if (succ != 0)
                tail_dist[op] = max(tail_dist[op], tail_dist[succ]);
        }
        tail_dist[op] += my_inst->time_idx[op];
        fo_online = max(fo_online, tail_dist[op]);
    }
    return ok;
}