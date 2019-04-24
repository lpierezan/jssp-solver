#include "global.hpp"
#include "instance.hpp"
#include "permutation.hpp"
#include <iostream>
#include <algorithm>
#include <utility>
#include <set>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

void instance::initialize(string filepath)
{
    // ===== Open instance file =================
    ifstream inst_file;
    inst_file.open(filepath, ios::in);
    inst_filepath = filepath;
    string tmp;
    getline(inst_file, tmp);
    inst_file >> n >> m;
    nop = n * m;

    getline(inst_file, tmp);
    getline(inst_file, tmp);

    //Init variables
    time_idx.resize(nop + 1, 0);
    mch_idx.resize(nop + 1, 0);
    js_idx.resize(nop + 1, 0);
    jp_idx.resize(nop + 1, 0);
    job_idx.resize(nop + 1, 0);

    time.resize(n + 1);
    machine.resize(n + 1);
    jm_to_pos.resize(n + 1);
    for (int i = 1; i <= n; i++)
    {
        time[i].resize(m + 1);
        machine[i].resize(m + 1);
        jm_to_pos[i].resize(m + 1);
    }

    // Operations Times
    int i, j;
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= m; j++)
        {
            inst_file >> time[i][j];
            time_idx[idx(i, j)] = time[i][j];
            job_idx[idx(i, j)] = i;
        }
    }

    inst_file >> tmp;

    //Operations Machines
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= m; j++)
        {
            inst_file >> machine[i][j];
            mch_idx[idx(i, j)] = machine[i][j];

            jm_to_pos[i][machine[i][j]] = j;
        }
    }

    inst_file.close();
    //===========

    build_job_neighborhood();
}

// ============ Algorithms ==============================

void instance::build_job_neighborhood()
{
    js_idx.resize(nop + 1, 0);
    jp_idx.resize(nop + 1, 0);

    for (int ji = 1; ji <= n; ji++)
    {
        for (int pi = 1; pi <= m; pi++)
        {
            int my_idx = idx(ji, pi);

            if (pi - 1 >= 1)
            {
                jp_idx[my_idx] = idx(ji, pi - 1);
            }

            if (pi + 1 <= m)
            {
                js_idx[my_idx] = idx(ji, pi + 1);
            }
        }
    }
}

void instance::build_mch_neighborhood(vvi &mch_perm, vector<int> &ms, vector<int> &mp)
{
    ms.resize(nop + 1, 0);
    mp.resize(nop + 1, 0);

    int i, mi, ji, jj, pos1, pos2, my_idx;
    for (mi = 0; mi < m; mi++)
    {

        for (i = 0; i < n; i++)
        {
            ji = mch_perm[mi][i];
            pos1 = jm_to_pos[ji][mi + 1];
            my_idx = idx(ji, pos1);
            //operation <ji,pos>
            mp[my_idx] = ms[my_idx] = 0;

            //machine pred.
            if (i - 1 >= 0)
            {
                jj = mch_perm[mi][i - 1];
                pos2 = jm_to_pos[jj][mi + 1];
                mp[my_idx] = idx(jj, pos2);
            }

            //machine succ.
            if (i + 1 < n)
            {
                jj = mch_perm[mi][i + 1];
                pos2 = jm_to_pos[jj][mi + 1];
                ms[my_idx] = idx(jj, pos2);
            }

            if (!check_mch_ng(my_idx, ms, mp))
            {
                cout << "Erro"
                        << ms[my_idx] << " "
                        << mp[my_idx] << " "
                        << pii_to_str(pair_job_mch(my_idx)) << " "
                        << pii_to_str(pair_job_mch(ms[my_idx])) << " "
                        << pii_to_str(pair_job_mch(mp[my_idx])) << " "
                        << endl;
                exit(1);
            }
        }
    }
}

bool instance::dfs(operation op, vector<vector<char>> &mark, vi &ms, vector<operation> &tsv)
{

    mark[op.ji][op.pi] = 1;

    operation ng;
    bool ok;

    //job succ.
    if (op.pi < m)
    {
        ng = operation(op.ji, op.pi + 1);
        switch (mark[ng.ji][ng.pi])
        {
        case 0:
            ok = dfs(ng, mark, ms, tsv);
            if (!ok)
                return false;
            break;
        case 1:
            cout << "cycle found" << endl;
            return false;
        case 2:
        default:
            break;
        }
    }

    //machine succ.
    int ms_idx = ms[idx(op)];
    if (ms_idx != 0)
    {
        auto op_pair = pair_job_pos(ms_idx);
        ng = operation(op_pair.first, op_pair.second);

        switch (mark[ng.ji][ng.pi])
        {
        case 0:
            ok = dfs(ng, mark, ms, tsv);
            if (!ok)
                return false;
            break;
        case 1:
            cout << "cycle found" << endl;
            return false;
        case 2:
        default:
            break;
        }
    }

    mark[op.ji][op.pi] = 2;
    tsv.push_back(op);
    return true;
}

bool instance::compute_fo(vvi mch_perm, int &fo)
{
    vi ms, mp;
    //======= Build DG Neighborhood
    build_mch_neighborhood(mch_perm, ms, mp);

    //======= Do topological sort
    vector<vector<char>> mark;
    vector<operation> tsv;
    mark.resize(n + 1);
    for (int i = 1; i <= n; i++)
        mark[i].resize(m + 1, 0);

    for (int ji = 1; ji <= n; ji++)
    {
        operation op(ji, 1);
        if (!mark[op.ji][op.pi])
        {
            auto ok = dfs(op, mark, ms, tsv);
            if (!ok)
                return false;
        }
    }

    if (tsv.size() != n * m)
    {
        cout << "not connected graph" << endl;
        return false;
    }
    reverse(ALL(tsv));

    //===== Compute longest path
    int makespan = 0;
    vi dist(n * m + 1, 0);
    for (int i = 0; i < tsv.size(); i++)
    {
        auto op = tsv[i];
        int my_idx = idx(op);

        dist[my_idx] = 0;
        int idx_prev;
        //job pred
        if (op.pi - 1 >= 1)
        {
            idx_prev = idx(op.ji, op.pi - 1);
            dist[my_idx] = max(dist[my_idx], dist[idx_prev]);
        }
        //mch pred
        if (mp[my_idx] != 0)
        {
            idx_prev = mp[my_idx];
            dist[my_idx] = max(dist[my_idx], dist[idx_prev]);
        }

        dist[my_idx] += time[op.ji][op.pi];
        makespan = max(makespan, dist[my_idx]);
    }

    fo = makespan;
    return true;
}

//Build active scheduling given operation order. O(n*m*log(n))
void instance::build_active_scheduling(vi &op_perm, vvi &mch_perm)
{
    //===================
    //Geenrating Active Schedule
    vector<set<pair<int, int>>> mch_sets(m);
    vi begin_time(nop + 1, -1);
    int op;
    for (int i = 0; i < op_perm.size(); i++)
    {
        op = op_perm[i];
        int jpred = jp_idx[op];
        //op must begin after pred_cut
        int pred_cut = 0;
        if (jpred != 0)
        {
            pred_cut = begin_time[jpred] + time_idx[jpred];
        }

        //loging for machine gaps
        int mch = mch_idx[op];
        int op_time = time_idx[op];
        mch--;

        if (mch_sets[mch].size() == 0)
        {
            //first job in machine
            begin_time[op] = pred_cut;
            mch_sets[mch].insert(make_pair(pred_cut, op));
        }
        else
        {
            //looking for mch cut
            auto mch_it = mch_sets[mch].upper_bound(make_pair(pred_cut, nop + 5)); //O(log(n))
            int mch_cut = pred_cut;
            if (mch_it != mch_sets[mch].begin())
            {
                auto mch_prev = --mch_it;
                mch_it++;
                mch_cut = mch_prev->first + time_idx[mch_prev->second];
            }

            //trying placements
            int try_time = max(mch_cut, pred_cut);
            while (mch_it != mch_sets[mch].end() && (try_time + op_time > mch_it->first))
            {
                try_time = mch_it->first + time_idx[mch_it->second];
                mch_it++;
            }

            begin_time[op] = try_time;
            mch_sets[mch].insert(make_pair(try_time, op));
        }
    }

    // =========================
    // Generating Machine Permutations
    mch_perm.clear();
    mch_perm.resize(m);
    for (int mi = 0; mi < m; mi++)
    {
        mch_perm[mi].resize(n);
        int pi = 0;
        for (auto op_it : mch_sets[mi])
        {
            mch_perm[mi][pi] = pair_job_pos(op_it.second).first;
            pi++;
        }
    }
}

// ================== Utilities =============================
void instance::print()
{
    cout << "n = " << n << " m = " << m << endl;
    cout << "Machine(Time)" << endl;

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            cout << setw(3) << machine[i][j] << "(" << setw(3) << time[i][j] << ") ";
        }
        cout << endl;
    }
}

string instance::mch_perm_to_str(const vvi &mch_perm)
{
    stringstream ss;
    ss << "Machine Permutations (of jobs):" << endl;
    for (int i = 0; i < m; i++)
    {
        ss << "machine " << (i + 1) << " job order: ";

        for (int j = 0; j < n; j++)
        {
            ss << (mch_perm[i][j]) << " ";
        }
        ss << endl;
    }
    ss << "==============" << endl;

    return ss.str();
}

bool instance::hard_check(vvi const &machine_perm, int cmp_fo)
{
    if (machine_perm.size() != m)
        return false;
    for (int mi = 0; mi < m; mi++)
    {
        if (machine_perm[mi].size() != n)
            return false;

        //permutation of [1...n]
        if (!permutation::check(machine_perm[mi], 1, n))
        {
            cout << "permutation test fail." << endl;
            return false;
        }
    }

    int fo;
    auto ok = compute_fo(machine_perm, fo);

    if (!ok)
    {
        cout << "Could not compute fo. Possible Cycle." << endl;
        return false;
    }
    if (cmp_fo != I32_INF && cmp_fo != fo)
    {
        cout << "check fo = " << fo << " != " << cmp_fo << endl;
        return false;
    }
    return true;
}

void instance::job_perm_to_op_perm(vi &job_perm)
{
    vi pi_cout(n + 1, 0);
    for (int i = 0; i < job_perm.size(); i++)
    {
        int ji = job_perm[i];
        int pi = pi_cout[ji] + 1;
        pi_cout[ji]++;
        job_perm[i] = idx(ji, pi);
    }
}

void instance::restart_pi_gen(int seed)
{
    my_rand_gen = std::default_random_engine(seed);
}

vi instance::generate_random_pi(int min, int max)
{
    return permutation::generate_rand(min, max, my_rand_gen);
}

string instance::op_vector_to_str(vi &op_vec, string sep)
{
    string ret;
    for (int i = 0; i < op_vec.size(); i++)
    {
        ret += pii_to_str(pair_job_pos(op_vec[i]));
        if (i < op_vec.size() - 1)
            ret += sep;
    }
    return ret;
}

bool instance::check_mch_ng(int v, vi &ms, vi &mp)
{
    auto pv = pair_job_mch(v);
    int prox = ms[v];
    auto pp = pair_job_mch(prox);
    int ant = mp[v];
    auto pa = pair_job_mch(ant);

    if (prox != 0)
    {
        if (pv.first == pp.first || pv.second != pp.second)
            return false;
    }

    if (ant != 0)
    {
        if (pv.first == pa.first || pv.second != pa.second)
            return false;
    }
    return true;
}