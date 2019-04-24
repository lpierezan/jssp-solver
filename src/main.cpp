#include <chrono>
#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <list>
#include <fstream>
#include <string>
#include <algorithm>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <random>
#include <functional>
#include <tuple>
#include <typeinfo>
#include <iomanip>
#include <sstream>
#include "global.hpp"
#include "permutation.hpp"
#include "time_control.hpp"
#include "stats_stream.hpp"

using namespace std;

// ============== Instance and Solution representation ===============

struct operation
{
public:
    int ji, pi;
    operation(){};
    operation(int _ji, int _pi) : ji(_ji), pi(_pi){};
};

class instance
{
    string inst_filepath;
    std::default_random_engine my_rand_gen = std::default_random_engine(0);

public:
    int n, m, nop;
    vvi time, machine; //[job idx][pos idx]
    vvi jm_to_pos;

    vi time_idx, mch_idx, js_idx, jp_idx, job_idx;

    // ==================== Initialization ==================
    void initialize(string filepath)
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

    void build_job_neighborhood()
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

    void build_mch_neighborhood(vvi &mch_perm, vector<int> &ms, vector<int> &mp)
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

    bool dfs(operation op, vector<vector<char>> &mark, vi &ms, vector<operation> &tsv)
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

    bool compute_fo(vvi mch_perm, int &fo)
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
    void build_active_scheduling(vi &op_perm, vvi &mch_perm)
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
    void print()
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

    string mch_perm_to_str(const vvi &mch_perm)
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

    inline pair<int, int> pair_job_pos(int idx)
    {
        int pi = idx % m;
        if (pi == 0)
            pi = m;
        int ji = ((idx - pi) / m) + 1;
        return make_pair(ji, pi);
    }

    inline pair<int, int> pair_job_mch(int idx)
    {
        return make_pair(job_idx[idx], mch_idx[idx]);
    }

    inline int idx(int ji, int pi)
    {
        return (ji - 1) * m + pi;
    }

    inline int idx(operation op)
    {
        return idx(op.ji, op.pi);
    }

    bool hard_check(vvi const &machine_perm, int cmp_fo = I32_INF)
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

    void job_perm_to_op_perm(vi &job_perm)
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

    void restart_pi_gen(int seed = 0)
    {
        my_rand_gen = std::default_random_engine(seed);
    }

    vector<int> generate_random_pi(int min, int max)
    {
        return permutation::generate_rand(min, max, my_rand_gen);
    }

    string op_vector_to_str(vi &op_vec, string sep = " ")
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

    bool check_mch_ng(int v, vi &ms, vi &mp)
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
};

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

    solution()
    {
        valid = false;
        my_inst = 0;
    };
    solution(instance &inst)
    {
        associate(inst);
    }

    void associate(instance &inst)
    {
        my_inst = &inst;
        n = my_inst->n;
        m = my_inst->m;
        nop = my_inst->nop;
        valid = false;
    }

    bool config(vvi &mch_perm)
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

    bool config(vi &op_perm)
    {
        my_mch_perm.clear();
        my_inst->build_active_scheduling(op_perm, my_mch_perm);
        return config(my_mch_perm);
    }

    bool repair()
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

    int get_fo()
    {
        return fo_online;
    }

    // ============ Utilities ======================
    int inline mch_idx(int op) { return my_inst->mch_idx[op]; }

    pair<int, int> get_succ(int op)
    {
        return make_pair(my_inst->js_idx[op], ms_idx[op]);
    }

    pair<int, int> get_pred(int op)
    {
        return make_pair(my_inst->jp_idx[op], mp_idx[op]);
    }

    int inline js(int idx) { return my_inst->js_idx[idx]; }

    int inline jp(int idx) { return my_inst->jp_idx[idx]; }

    int inline idx(int ji, int pi)
    {
        return my_inst->idx(ji, pi);
    }

    bool inline is_critical(int op)
    {
        //0 -> op + op -> end
        return ((tail_dist[op] + head_dist[op] - my_inst->time_idx[op]) == fo_online);
    }

    bool inline is_block_start(int op)
    {
        return (is_critical(op) && (mp_idx[op] == 0 || !is_critical(mp_idx[op])));
    }

    void print(bool print_mach_perm = false)
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

private:
    //Dfs used for topological sort
    bool ts_dfs(int op_idx)
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
    vector<char> mark;
    bool build_ts()
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
    bool build_dists()
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
};

class solution_generator
{
    instance *my_inst;
    std::default_random_engine my_rand_gen = std::default_random_engine(0);

public:
    solution_generator() {}
    solution_generator(instance &inst, std::default_random_engine _rand_gen = std::default_random_engine(0))
    {
        my_inst = &inst;
        my_rand_gen = _rand_gen;
    }

    void reset_rand(int seed = 0)
    {
        my_rand_gen = std::default_random_engine(seed);
    }

    solution generate()
    {
        solution s(*my_inst);

        // == Permutation of all operations (Using repeted job_index permutation) =====
        vi job_rep_perm;
        for (int ji = 1; ji <= my_inst->n; ji++)
        {
            for (int pi = 1; pi <= my_inst->m; pi++)
            {
                job_rep_perm.push_back(ji);
            }
        }
        //Shuffle
        shuffle(ALL(job_rep_perm), my_rand_gen);

        //From job index => operation index
        my_inst->job_perm_to_op_perm(job_rep_perm);

        // ======= Config Solution ===========
        auto ok = s.config(job_rep_perm);
        if (!ok)
        {
            cout << "Error - Invalid solution generated. (not supposed to happen)" << endl;
            exit(0);
        }

        return s;
    }
};
// =================== Local Search ========================================

class ls_engine
{
public:
    // ======= Meta data  ============
    struct move
    {
        int op1, op2;
        bool fw;

        string to_string(instance *inst)
        {
            return pii_to_str(inst->pair_job_pos(op1)) + " " + (fw ? "fw" : "bw") + " " +
                   pii_to_str(inst->pair_job_pos(op2));
        }

        bool check(solution &s)
        {
            auto p1 = s.my_inst->pair_job_mch(op1);
            auto p2 = s.my_inst->pair_job_mch(op2);

            if (p1.first == p2.first)
                return false;
            if (p1.second != p2.second)
                return false;

            return true;
        }
    };

    class move_mem_data
    {
        vi block_idxs;
        static std::hash<int> hsh_int;

    public:
        bool fw;
        size_t seed = 0;

        void init(bool _fw)
        {
            fw = _fw;
            int x_fw = fw ? 2 : 1;
            seed ^= hsh_int(x_fw) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

        struct eqp_cmp
        {
            bool operator()(const move_mem_data &move_mem1, const move_mem_data &move_mem2) const
            {
                if (move_mem2.fw != move_mem1.fw)
                    return false;
                if (move_mem1.block_idxs.size() != move_mem2.block_idxs.size())
                    return false;
                for (int i = 0; i < move_mem1.block_idxs.size(); i++)
                    if (move_mem1.block_idxs[i] != move_mem2.block_idxs[i])
                        return false;
                return true;
            }
        };

        void add_value(int x)
        {
            //compund hash similar to boost implementation.
            seed ^= hsh_int(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            block_idxs.push_back(x);
        }

        struct move_hash
        {
            size_t operator()(const move_mem_data &move) const
            {
                return move.seed;
            }
        };

        string to_string()
        {
            string ret = (fw ? "fw " : "bw ");
            for (int x : block_idxs)
            {
                ret += std::to_string(x) + " ";
            }
            return ret;
        }
    };

    // ====== Memory data =========
    struct ls_stats
    {
        stats_stream<long> one_round_time_stats = stats_stream<long>("round time [ms]");
        long total_time = 0, n_ls_moves = 0;
        vector<move> moves_memory;
        long n_moves_no_local_imp = 0;
        long n_moves_no_global_imp = 0;

        vvi best_sol_found; //machine perm representation
        int best_fo_found = I32_INF;
        bool stop_because_no_move = false;

        unordered_set<move_mem_data, move_mem_data::move_hash, move_mem_data::eqp_cmp> tabu_set;
        deque<move_mem_data> tabu_list;

        void reset(solution &s)
        {
            total_time = n_ls_moves = n_moves_no_local_imp = 0;
            n_moves_no_global_imp = 0;
            moves_memory.clear();
            one_round_time_stats = stats_stream<long>("round time [ms]");
            best_sol_found = s.my_mch_perm; //copy
            tabu_set.clear();
            tabu_list.clear();
            best_fo_found = I32_INF;
            if (s.valid)
                best_fo_found = s.fo_online;
            stop_because_no_move = false;
        }

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

        void compute(solution &s)
        {
            tabu_list_size = s.n * s.m * (tabu_list_factor >= 0 ? tabu_list_factor : 1);
        }

        ls_parameters() {}

    } params;

private:
public:
    // ============
    ls_engine()
    {
    }

    //Parameters Input
    void config(ls_parameters _params = ls_parameters())
    {
        params = _params;
    }

    //String signature
    string signature(char sep = ' ')
    {
        //TODO
        string ret = "";
        return ret;
    }

    void run_ls(solution &initial_sol)
    {
        //Setting up memory and stats variables
        stats.reset(initial_sol);
        //Initializing aux data used by the ls algorithm
        n = initial_sol.n;
        m = initial_sol.m;
        nop = initial_sol.nop;
        s = &initial_sol;
        mark.resize(nop + 1, false);
        aux_buf.resize(nop + 1);

        //Setting up ls parameters
        params.compute(initial_sol);

        //run ls
        ls_start(initial_sol);

        //updating the given initial solution to the best
        initial_sol.config(stats.best_sol_found);
    }

    int best_fo_found()
    {
        return stats.best_fo_found;
    }

    //returns the job permutation in each machine
    vvi best_sol_found()
    {
        //not neccessarly during the local search
        return stats.best_sol_found;
    }

private:
    // ======== Locas Search Algorithms =============
    int n, m, nop;
    solution *s;
    vi hb_dist, tb_dist; //block distances
    vector<bool> mark;   //aux mark vector
    vector<int> aux_buf;
    vi explored_blocks;

    bool tabu_acc(move_mem_data &move_tabu_data, int eval)
    {
        //Accept improve moves or not tabu moves.
        bool tabu_acc = false;
        if (eval < best_fo_found())
            return true;

        //query tabu
        auto it = stats.tabu_set.find(move_tabu_data);
        tabu_acc = (it == stats.tabu_set.end());

        if (params.debug_mode)
        {
            cout << "Tabu " << (tabu_acc ? "Acc" : "Rej") << " move: "
                 << move_tabu_data.to_string()
                 << endl;
        }

        return tabu_acc;
    }

    void tabu_add(move_mem_data &move_tabu_data)
    {
        auto &tabu_list = stats.tabu_list;
        stats.tabu_set.insert(move_tabu_data);
        tabu_list.push_back(move_tabu_data);

        if (params.debug_mode)
        {
            cout << "Inserting Tabu Move: " + move_tabu_data.to_string() << endl;
        }

        /* === Tabu list management: ===
		If the list is full => Remove the oldest move
		*/
        while (tabu_list.size() > params.tabu_list_size)
        {
            auto move_remove = tabu_list.front();
            tabu_list.pop_front();
            stats.tabu_set.erase(move_remove);
            if (params.debug_mode)
            {
                cout << "Removing Tabu Move: " + tabu_list.front().to_string() << endl;
            }
        }
    }

    void create_tabu_move_data(move &move, vi &ms, vi &mp,
                               vector<move_mem_data> &v_tabu_move_data)
    {
        /* Create a tabu move data related to the move u fw v.
		The algorithm is used in a symmetric manner to compute u bw v*/

        /* <u,x, .... ,v..>  => <x,...,v,u,...>
		u fw v => tabu u bw x (Correct but may not exist) (*)
		u fw v => tabu x fw u (Strong Tabu Rule)
		u fw v in <u,v,...> => tabu v fw u

		Break of symmetry:
		(*) when the move is u bw v and move = <u b> it is named v fw u.
		u bw v => tabu x bw v => tabu x fw v

		*/

        int u = move.op1;
        int v = move.op2;
        int x = ms[u], cur;
        move_mem_data tabu1, tabu2;
        tabu1.init(!move.fw);
        tabu1.add_value(u);
        tabu1.add_value(v);
        cur = v;
        while (cur != x)
        {
            cur = mp[cur];
            tabu1.add_value(cur);
        }
        v_tabu_move_data.push_back(tabu1);

        if (v == ms[u])
        { //Case (1) move = <u,v...>
            //v fw u
            tabu2.init(move.fw);
            tabu2.add_value(v);
            tabu2.add_value(u);
            v_tabu_move_data.push_back(tabu2);
        }
        else
        { //Case (2) move = <u,x...,v>
            if (params.strong_tabu_rule)
            {
                //x fw u
                tabu2.init(move.fw);
                tabu2.add_value(x);
                cur = x;
                while (cur != v)
                {
                    cur = ms[cur];
                    tabu2.add_value(cur);
                }
                tabu2.add_value(u);
                v_tabu_move_data.push_back(tabu2);
            }
        }
    }

    void tabu_list_size_update()
    {
        if (!params.change_tabu_list_size)
            return;

        int change_max = params.tabu_list_percent_change; //%
        if ((stats.n_ls_moves % (params.max_moves_no_global_improve / 10)) == 0)
        {
            int x = (rand_gen() % (2U * change_max));
            x -= change_max;

            double percent_change = 1.0 + x / 100.0;
            params.tabu_list_size *= percent_change;
            params.tabu_list_size = max(params.tabu_list_size, 5);
        }
    }

    bool test_acyclic_move(int u, int v, bool is_fw)
    {
        //Reference Zhang2007
        bool is_acy = false;
        if (is_fw)
        { // u in front of v
            //Ensure that does not exist path js[u] -> v.
            int u_suc = s->js(u);
            is_acy = s->tail_dist[v] >= s->tail_dist[u_suc];
        }
        else
        { // u before v
            //Ensure that does not exist path v -> jp[u].
            int u_pred = s->jp(u);
            is_acy = s->head_dist[v] >= s->head_dist[u_pred];
        }

        return is_acy;
    }

    void explore_block_moves(vi &block_idxs,
                             vi &js, vi &jp, vi &ms, vi &mp, vi &head_dist, vi &tail_dist, vi &time,
                             int &best_eval, move &best_move,
                             bool fw)
    {
        /* The code is based on the move  u -> v.
		But, this method is called also for move v <- u with the symmetric
		parameters.
		*/

        int block_size = block_idxs.size();

        int u = block_idxs[0];
        int pu = time[u];
        int u_js = js[u];
        int u_jp = jp[u];

        move_mem_data move_tabu_data;
        move_tabu_data.init(fw);
        move_tabu_data.add_value(u);

        for (int i = 1; i < block_size; i++)
        {
            int v = block_idxs[i];
            move_tabu_data.add_value(v);
            //===========================
            // Longest distance approx from 0 -> v (after remove of vertex u from the beginning of block)
            int vpred = jp[v], vsuc = js[v], pv = time[v];
            if (i == 1)
            {
                hb_dist[i] = head_dist[mp[u]];
            }
            else
            {
                hb_dist[i] = hb_dist[i - 1];
            }
            hb_dist[i] = max(hb_dist[i], head_dist[vpred]); //h(0) = 0
            hb_dist[i] += pv;
            //==========================

            if (params.debug_mode)
            {
                cout << "Test move " << pii_to_str(s->my_inst->pair_job_pos(u)) << " "
                     << (fw ? "fw" : "bw") << " " << pii_to_str(s->my_inst->pair_job_pos(v)) << endl;
            }

            //Test if move u->v is acyclic. If its not, ignore move.
            bool acly_move = test_acyclic_move(u, v, fw);
            if (!acly_move)
            {
                //If a move u->v is not acyclic the next move isn�t too.
                if (params.debug_mode)
                    cout << "Move can make cycle." << endl;

                break;
            }
            if (params.debug_mode)
            {
                cout << "Move not makes cycles." << endl;
            }

            //=========== EVALUATION OF MOVE U -> V ================

            // === approx distance from 0 -> u
            hb_dist[0] = hb_dist[i];
            hb_dist[0] = max(hb_dist[0], head_dist[u_jp]); //h(0) = 0
            hb_dist[0] += pu;
            // ===

            // === approx dist from u -> end
            tb_dist[0] = tail_dist[u_js];                   //t(0) = 0
            tb_dist[0] = max(tb_dist[0], tail_dist[ms[v]]); //t(0) = 0
            tb_dist[0] += pu;
            // ===

            //Compute Longest distance approx from v� [v� <= v]  (after insert of vertex u after v)
            // === approx dist from v -> end
            tb_dist[i] = tb_dist[0];
            tb_dist[i] = max(tb_dist[i], tail_dist[vsuc]); //t(0) = 0
            tb_dist[i] += pv;
            // ===

            int j, tmp;
            for (j = i - 1; j >= 1; j--)
            {
                int a = block_idxs[j];
                int asuc = js[a];

                tmp = tb_dist[j];
                tb_dist[j] = tb_dist[j + 1];
                tb_dist[j] = max(tb_dist[j], tail_dist[asuc]); //t(0) = 0
                tb_dist[j] += time[a];
            }

            // ==== evaluate move u -> v ====
            int eval = 0;
            for (j = 0; j <= i; j++)
            {
                int p_a = time[block_idxs[j]];
                eval = max(eval, hb_dist[j] + tb_dist[j] - p_a);
            }

            //=====

            if (params.debug_mode)
            {
                cout << "Move evaluation = " << eval << endl;
            }

            //Testing if evaluation is the new best
            bool break_even = params.break_even_moves;
            bool is_new_best = (best_eval == -1 || (eval < best_eval) || (break_even && (eval == best_eval) && s->is_critical(js[v])));

            if (is_new_best)
            {
                //=====================================
                //Test if move is tabu accepted
                bool tabu_accepted = false;
                tabu_accepted = tabu_acc(move_tabu_data, eval);
                //====================================
                if (tabu_accepted)
                {

                    best_move.op1 = u;
                    best_move.op2 = v;
                    best_move.fw = fw;
                    best_eval = eval;

                    /*
					if (!best_move.check(*s)) {
					cout << "Error invalid move found." << endl;
					exit(1);
					}*/

                    if (params.debug_mode)
                    {
                        cout << "New best." << endl;
                    }
                }
                else
                {
                    if (params.debug_mode)
                    {
                        cout << "Tabu Rejected." << endl;
                    }
                }
            }
        }
    }

    bool search_block(int block_begin, move &best_move, int &best_eval)
    {
        //return true to stop dfs
        if (block_begin == 0 || mark[block_begin])
        {
            return false;
        }

        mark[block_begin] = true;
        explored_blocks.push_back(block_begin);

        int block_end = block_begin, mch_succ;
        vi block_idxs;
        block_idxs.push_back(block_begin);

        mch_succ = s->ms_idx[block_end];
        while (mch_succ != 0 && s->is_critical(mch_succ))
        {
            block_end = mch_succ;
            mch_succ = s->ms_idx[block_end];
            block_idxs.push_back(block_end);
        }

        if (params.debug_mode)
        {
            cout << "critical block: ";
            cout << s->my_inst->op_vector_to_str(block_idxs, " -> ") << endl;
        }

        //initializing aux vector (block distance vectors)
        if (block_idxs.size() > hb_dist.size())
        {
            hb_dist.resize(block_idxs.size());
            tb_dist.resize(block_idxs.size());
        }

        if (block_idxs.size() > 1)
        {

            //explore moves u-> foward
            explore_block_moves(block_idxs, s->my_inst->js_idx, s->my_inst->jp_idx, s->ms_idx, s->mp_idx,
                                s->head_dist, s->tail_dist, s->my_inst->time_idx, best_eval, best_move, true);

            if (block_idxs.size() > 2)
            {
                //explore moves backward <- v (the same algorithm but symmetric)
                reverse(ALL(block_idxs));
                explore_block_moves(block_idxs, s->my_inst->jp_idx, s->my_inst->js_idx, s->mp_idx, s->ms_idx,
                                    s->tail_dist, s->head_dist, s->my_inst->time_idx, best_eval, best_move, false);
                reverse(ALL(block_idxs));
            }
        }

        if (params.debug_mode)
        {
            cout << "Best move so far: " << best_move.to_string(s->my_inst)
                 << " Best eval so far = " << best_eval << endl;
        }

        if (best_eval != -1)
        {
            if (params.stop_first_improve && best_eval < s->fo_online)
            {
                return true; // stop search
            }
        }

        //dfs to next blocks on the same path (job successor move)
        bool stop_search = false;
        int js_block_end = s->my_inst->js_idx[block_end];
        if (!params.explore_all_blocks)
            stop_search = search_block(js_block_end, best_move, best_eval);
        else if (s->is_block_start(js_block_end))
            stop_search = search_block(js_block_end, best_move, best_eval);

        if (stop_search)
            return true;

        //dfs to parallel blocks
        if (params.explore_all_blocks)
        {
            for (int i = 0; i < block_idxs.size(); i++)
            {
                int job_succ = s->my_inst->js_idx[block_idxs[i]];
                if (s->is_block_start(job_succ))
                {
                    stop_search = search_block(job_succ, best_move, best_eval);
                }
                if (stop_search)
                    return true;
            }
        }
        return stop_search;
    }

    bool search_move(move &selected_move)
    {
        bool found_move = false;

        /*
		Search for moves:
		1 - Scan critical blocks of size > 1
		2 - Look for fw and bw interchanges of N7 (Balas1998 + Zhang2007)
		3 - Evaluate using approx. method (Balas 1998)
		4 - Test if its Tabu Accepted or Improve move

		Returns the best approx. evaluated move that is improve or tabu accepted.
		*/

        //scan blocks
        if (params.debug_mode)
        {
            cout << "Beginning Block Scan" << endl;
        }

        vector<pair<move, int>> possible_moves;
        int i, ji, pi, mi, op, block_begin = 0;
        vi block_start_options, block_order;

        //find beginning of first blocks
        for (mi = 1; mi <= m; mi++)
        {
            ji = s->my_mch_perm[mi - 1][0];
            pi = s->my_inst->jm_to_pos[ji][mi];

            op = s->idx(ji, pi);
            if (s->is_block_start(op))
            {
                block_start_options.push_back(op);
                block_order.push_back(block_start_options.size() - 1);
                if (!params.choose_rand_first_block && !params.explore_all_blocks)
                    break;
            }
        }

        if (params.choose_rand_first_block)
            shuffle(ALL(block_order), rand_gen);

        //perform dfs block exploration
        move best_move;
        int best_eval = -1;

        explored_blocks.clear();
        for (i = 0; i < block_order.size(); i++)
        {
            block_begin = block_start_options[block_order[i]];
            search_block(block_begin, best_move, best_eval);
            if (!params.explore_all_blocks)
                break;
        }

        //clean dfs of blocks mark
        for (i = 0; i < explored_blocks.size(); i++)
            mark[explored_blocks[i]] = false;

        if (best_eval != -1)
        {
            selected_move = best_move;

            if (params.not_down_moves && best_eval > s->fo_online)
            {
                found_move = false;
                if (params.debug_mode)
                {
                    cout << "Down Move - Critical Path Move Rejected: " << selected_move.to_string(s->my_inst) << " eval = " << best_eval << endl;
                }
            }
            else
            {
                found_move = true;

                if (params.debug_mode)
                {
                    cout << "Critical Path Move Selected: " << selected_move.to_string(s->my_inst) << " eval = " << best_eval << endl;
                }
            }
        }
        else
        {
            if (params.debug_mode)
            {
                cout << "No move found in critical path." << endl;
            }
            found_move = false;
        }

        return found_move;
    }

    void update_tsv(int u, int v, bool fw,
                    int &u_init_pos, int &v_init_pos)
    {
        //=== Updating topological sorted vector
        //Adapting Nowicki2005 - Head and Tail speed up

        if (params.debug_mode)
        {
            cout << "Updating Top. Sort Vector" << endl;
            cout << "tsv: " << s->my_inst->op_vector_to_str(s->tsv) << endl;
        }

        int inc;
        vi *js, *ms;
        if (fw)
        {
            inc = 1;
            js = &(s->my_inst->js_idx);
            ms = &(s->ms_idx);
        }
        else
        {
            inc = -1;
            js = &(s->my_inst->jp_idx);
            ms = &(s->mp_idx);
        }

        int i, j;
        vi &tsv = (s->tsv), &tsv_inv = (s->tsv_inv);

        int u_pos = tsv_inv[u], v_pos = tsv_inv[v];
        u_init_pos = u_pos;
        v_init_pos = v_pos;
        int nv = 0;
        int u_seq = v_pos;

        //including vertices x such that has path to v
        mark[v] = true;
        aux_buf[nv++] = v;
        for (i = v_pos - inc; i != (u_pos - inc); i -= inc)
        {
            int x = tsv[i];
            if (x != u && (mark[js->operator[](x)] || mark[ms->operator[](x)]))
            {
                mark[x] = true;
                aux_buf[nv++] = x;
            }
            else
            {
                tsv[u_seq] = x;
                u_seq -= inc;
            }
        }
        for (i = nv - 1, j = u_pos; i >= 0; i--, j += inc)
        {
            tsv[j] = aux_buf[i];
            mark[aux_buf[i]] = false;
        }
        //tsv: <x1,x2...,v> <u,....>	, for fw moves

        //update tsv_inv
        for (i = u_pos; i != (v_pos + inc); i += inc)
            tsv_inv[tsv[i]] = i;

        if (params.debug_mode)
        {
            cout << "Updaed Top. Sort Vector" << endl;
            cout << "tsv: " << s->my_inst->op_vector_to_str(s->tsv) << endl;
        }
    }

    void apply_move(move &move_apply)
    {
        /*
		To apply a move we need to update the solution structure:
		- mp and ms ( O(1) )
		- topological sort vector (O(nop))*
		- head and tail dists ( O(nop) - can improve? )
		- fo (O(1) given head and tail dists)
		- machine permutation ( easy O(m) - more complex O(|block|)
		*Adapting Nowicki2005 - Head and Tail speed up

		AND
		Add the move in tabu.
		*/

        if (params.debug_mode)
        {
            cout << "Applying move." << endl;
        }
        vi &tsv = s->tsv, &time = s->my_inst->time_idx,
           &js = s->my_inst->js_idx, &jp = s->my_inst->jp_idx,
           &ms = s->ms_idx, &mp = s->mp_idx,
           &hd = s->head_dist, &td = s->tail_dist;

        int u = move_apply.op1;
        int v = move_apply.op2;
        bool fw = move_apply.fw;
        int inc = fw ? 1 : -1;

        // === Including move in tabu
        vector<move_mem_data> v_tabu_move_data;
        create_tabu_move_data(move_apply, fw ? ms : mp, fw ? mp : ms, v_tabu_move_data);
        for (auto tabu_move_data : v_tabu_move_data)
            tabu_add(tabu_move_data);

        tabu_list_size_update();
        //====

        //=== Update topological sorted vector
        int u_pos, v_pos;
        update_tsv(u, v, fw, u_pos, v_pos);

        //=== Update mp,ms
        if (mp[u] != 0)
            ms[mp[u]] = ms[u];
        if (ms[u] != 0)
            mp[ms[u]] = mp[u];
        if (fw)
        { //u in front of v
            mp[u] = v;
            ms[u] = ms[v];
            ms[v] = u;
            if (ms[u] != 0)
                mp[ms[u]] = u;
        }
        else
        { //u before v
            ms[u] = v;
            mp[u] = mp[v];
            mp[v] = u;
            if (mp[u] != 0)
                ms[mp[u]] = u;
        }

        //=== Update distances
        int r_idx = max(u_pos, v_pos), l_idx = min(u_pos, v_pos);
        int i, x, &fo_online = (s->fo_online);

        //[l_idx,end] -> update head
        fo_online = 0; //update fo
        for (i = l_idx; i < tsv.size(); i++)
        {
            x = tsv[i];
            hd[x] = max(hd[jp[x]], hd[mp[x]]) + time[x];
            fo_online = max(fo_online, hd[x]);
        }

        //[0,r_idx] : update tail
        for (i = r_idx; i >= 0; i--)
        {
            x = tsv[i];
            td[x] = max(td[js[x]], td[ms[x]]) + time[x];
            fo_online = max(fo_online, hd[x]);
        }

        //=== Update machine perm
        int u_job = s->my_inst->job_idx[u],
            v_job = s->my_inst->job_idx[v],
            mch_row_idx = s->my_inst->mch_idx[u] - 1;
        vi &mch_prm = (s->my_mch_perm[mch_row_idx]);

        for (i = 0; mch_prm[i] != u_job; i++)
            ;
        while (mch_prm[i + inc] != v_job)
        {
            swap(mch_prm[i], mch_prm[i + inc]);
            i += inc;
        }
        swap(mch_prm[i], mch_prm[i + inc]);

        // ======================================

        if (params.debug_mode)
        {
            cout << "Move applied. Printing Solution: " << endl;
            s->print(true);

            if (!s->my_inst->hard_check(s->my_mch_perm, s->fo_online))
            {
                cout << "Hard Check Fail." << endl;
            }
        }
    }

    bool stop_criterion_test()
    {
        //Time limit
        long long dummy;
        if (tl_clock_check(dummy))
            return true;

        if (params.restrict_n_iterations && stats.n_ls_moves > params.max_n_moves)
            return true;

        if (params.restrict_max_no_global_improve && stats.n_moves_no_global_imp > params.max_moves_no_global_improve)
            return true;

        return false;
    }

    void keep_best_sol_test()
    {
        //Keep the best solution so far
        if (s->fo_online < stats.best_fo_found)
        {
            stats.best_fo_found = s->fo_online;
            stats.best_sol_found = s->my_mch_perm; //copy
            stats.n_moves_no_global_imp = 0;
        }
        else
        {
            stats.n_moves_no_global_imp++;
        }
    }

    void ls_start(solution &s)
    {

        auto start_ls_time = get_time::now();
        // ========== LS start ==============
        bool stop_ls = false;
        do
        {
            auto start_round_time = get_time::now();
            auto end_round_time = get_time::now();
            //======== round start =================

            //Search for good moves
            bool found_move = false;
            move selected_move;
            found_move = search_move(selected_move);

            if (!found_move)
            {
                //=========== stop counting round time
                end_round_time = get_time::now();
                //====================================
                stop_ls = true;
                stats.stop_because_no_move = true;
            }
            else
            {
                // ===== Apply Move =======
                auto old_fo = s.fo_online;
                apply_move(selected_move);
                keep_best_sol_test();
                auto new_fo = s.fo_online;

                if (new_fo < old_fo)
                {
                    //=========== Improve move =============
                    stats.n_moves_no_local_imp = 0;
                }
                else if (new_fo == old_fo)
                {
                    //=========== Side move ================
                    stats.n_moves_no_local_imp++;
                }
                else
                {
                    //=========== Unimprove move ================
                    stats.n_moves_no_local_imp++;
                }

                //=========== stop counting round time
                end_round_time = get_time::now();
                //====================================
            }
            stats.one_round_time_stats.add(to_time(end_round_time - start_round_time));

            //Stop
            if (stop_criterion_test())
            {
                stop_ls = true;
            }

            if (!stop_ls)
                stats.n_ls_moves++;
        } while (!stop_ls);

        // ===== end of ls ========
        auto end_ls_time = get_time::now();
        stats.total_time = to_time(end_ls_time - start_ls_time);
    }
};

// ================= Metaheuristics ====================

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

        mh_params()
        {
            // ============
            //Intensive TS - Random Solution Generated
            ls_params_indv_gen = ls_engine::ls_parameters();
            ls_params_indv_gen.max_moves_no_global_improve = 150;
            ls_params_indv_gen.tabu_list_factor = 1;
            // ============
            //Light TS - Path soution evaluation
            ls_params_path_eval = ls_engine::ls_parameters();
            ls_params_path_eval.max_moves_no_global_improve = 10;
            ls_params_path_eval.tabu_list_factor = 1;
            //============
            //Intensive TS - Offspring generated
            ls_params_offspring = ls_engine::ls_parameters();
            ls_params_offspring.max_moves_no_global_improve = 150;
            ls_params_offspring.tabu_list_factor = 1;
            //============
        }

        void init(instance &inst)
        {

            sol_generator = solution_generator(inst);
            if (curr_pop_size - n_elite - n_off_to_gen < 0)
            {
                cout << "Error (curr_pop_size - n_elite - n_off_to_gen < 0)";
                exit(1);
            }

            ls_params_indv_gen.max_moves_no_global_improve = inst.n * inst.m;
            ls_params_path_eval.max_moves_no_global_improve = inst.n * inst.m / break_path_blocks;
            ls_params_offspring.max_moves_no_global_improve = inst.n * inst.m;
        }
    } params;

    // ==== MH Memory/Stats Data
    struct mh_stats
    {
        //Solutions data
        int best_fo_found = I32_INF;
        solution best_sol_found; //only set at the end

        //Time measures
        long total_time = 0;
        chrono::time_point<chrono::steady_clock> start_time;
        stats_stream<long> it_time_stream;

        int curr_iteration = 0, n_iter_no_global_improve = 0;
        stats_stream<long long> parents_dist;
        int n_infeasible = 0;

        void init(instance &inst)
        {
            best_fo_found = I32_INF;
            best_sol_found = solution(inst);
            total_time = 0;
            start_time = get_time::now();
            it_time_stream.clear();
            curr_iteration = 0;
            n_iter_no_global_improve = 0;
            n_infeasible = 0;
        }

    } stats;

    // === Meta data ========
    class individual
    {
    public:
        vvi mch_perm;
        int fitness;

        individual(solution &s)
        {
            mch_perm = s.my_mch_perm;
            fitness = s.fo_online;
        }

        string to_string(instance &inst)
        {
            stringstream ret;
            ret << "Makespan = " << fitness << endl;
            ret << inst.mch_perm_to_str(mch_perm);
            return ret.str();
        }

        static individual make_invalid()
        {
            individual i;
            i.fitness = I32_INF;
            return i;
        }

    private:
        individual(){};
    };

    class population
    {
        vector<individual *> indv_pointers;
        int best_fit = I32_INF;
        bool sorted = false;

    public:
        ~population()
        {
            clear();
        }

        int best_fitness() { return best_fit; }

        int size()
        {
            return indv_pointers.size();
        }

        void clear()
        {
            for (auto p : indv_pointers)
                delete p;

            clear_no_delete();
        }

        void add(const individual &_indv)
        {
            sorted = false;
            individual *indv = new individual(_indv);
            indv_pointers.push_back(indv);
            best_fit = min(best_fit, indv->fitness);
        }

        //empty offspring
        void consume(population &off)
        {
            sorted = false;
            best_fit = min(best_fit, off.best_fit);
            indv_pointers.insert(indv_pointers.end(), ALL(off.indv_pointers));
            off.clear_no_delete();
        }

        void pop_back()
        {
            if (indv_pointers.size() == 0)
                return;
            auto p = indv_pointers[indv_pointers.size() - 1];
            delete p;
            indv_pointers.pop_back();
        }

        void sort()
        {

            std::sort(ALL(indv_pointers),
                      [](individual *p1, individual *p2) -> bool { return p1->fitness < p2->fitness; });

            sorted = true;
        }

        void keep_best(int n_best)
        {
            if (!sorted)
            {
                cout << "Not sorted error.";
                exit(1);
            }

            while (indv_pointers.size() > n_best)
            {
                pop_back();
            }
        }

        const individual *select_bests(int n_elite)
        {
            if (!sorted)
            {
                cout << "Not sorted.";
                exit(1);
            }
            n_elite = min(n_elite, (int)indv_pointers.size());
            if (n_elite == 0)
                return 0;
            int i = rand_gen() % n_elite;
            return indv_pointers[i];
        }

        const individual *select_worst(int n_worse)
        {
            if (!sorted)
            {
                cout << "Not sorted.";
                exit(1);
            }
            n_worse = min(n_worse, (int)indv_pointers.size());
            if (n_worse == 0)
                return 0;
            int i = rand_gen() % n_worse;
            i = (indv_pointers.size() - 1 - i); //i in [n-i-1, n-1]
            return (indv_pointers[i]);
        }

        string to_string(instance &inst)
        {
            stringstream ret;
            ret << "Bestfit = " << best_fit << endl;
            ret << "Size = " << size() << endl;
            ret << "Individuals: " << endl;
            for (auto p : indv_pointers)
                cout << p->to_string(inst) << endl;
        }

        individual get_best()
        {
            if (!sorted)
                sort();
            if (indv_pointers.size() == 0)
            {
                return individual::make_invalid();
            }
            else
            {
                return *(indv_pointers[0]);
            }
        }

    private:
        void clear_no_delete()
        {
            indv_pointers.clear();
            sorted = false;
            best_fit = I32_INF;
        }
    };

public:
    void config(mh_params _mh_params = mh_params())
    {
        params = _mh_params;
    }

    solution run(instance &inst_)
    {
        inst = &inst_;
        params.init(inst_);
        stats.init(inst_);
        //pre-allocation
        int n = inst->n;
        int m = inst->m;
        aux1.resize(m + 1);
        for (int i = 0; i <= m; i++)
            aux1[i].resize(n + 1, 0);
        //run
        run_mh();
        return stats.best_sol_found;
    }

    long get_best_fo() { return stats.best_fo_found; }

    solution get_best_sol()
    {
        return stats.best_sol_found;
    }

    string signature()
    { //TODO
        string ret = "";
        return ret;
    }

private:
    vvi aux1; //[m,n] pre-allocated

    individual gen_new_individual()
    {
        solution s = params.sol_generator.generate();
        //improve solution
        if (params.improve_after_generation)
        {
            params.ls.config(params.ls_params_indv_gen);
            params.ls.run_ls(s);
        }

        individual indv(s);
        if (params.log_mode)
        {
            cout << "New Indv Gen: " << endl;
            cout << indv.to_string(*inst) << endl;
        }

        return indv;
    }

    void gen_init_pop(population &pop)
    {
        /*
		Build a initial population.
		Individuals are generates and then local search is applyed.
		*/

        pop.clear();
        int i, n_pop = params.initial_pop_size;

        //configure ls algorithm for individual generation
        params.ls.config(params.ls_params_indv_gen);

        for (i = 0; i < n_pop; i++)
        {
            //generate new individual
            individual indv = gen_new_individual();
            //add to pop
            pop.add(indv);
        }

        //sorts the population (used to track elite solutions)
        pop.sort();
    }

    void update_best_fitness(population &pop)
    {
        auto pop_best = pop.best_fitness();
        if (pop_best < stats.best_fo_found)
        {
            //Improve makespan
            stats.n_iter_no_global_improve = 0;
        }
        else if (pop_best == stats.best_fo_found)
        {
            //Side move
            stats.n_iter_no_global_improve++;
        }
        else
        {
            cout << "Error. Makespan decreased.";
            exit(1);
        }

        stats.best_fo_found = pop_best;
    }

    bool stop_mh(population &pop)
    {

        //max iterations
        if (params.stop_max_it && stats.curr_iteration > params.max_it_limit)
            return true;
        //max no improve iterations
        if (params.stop_max_no_imp_it && stats.n_iter_no_global_improve > params.max_no_imp_it_limit)
            return true;
        //time limit
        long curr_time_sec = to_time(get_time::now() - stats.start_time) / 1000.0;
        if (params.stop_time_limit && curr_time_sec > params.time_limit_sec)
            return true;

        return false;
    }

    //must return <elite,normal>
    pair<const individual *, const individual *> select_parents(population &pop)
    {

        /*
		Methos based on
		A BIASED RANDOM-KEY GENETIC ALGORITHM FOR JOB-SHOP SCHEDULING
		JOS� FERNANDO GON�ALVES AND MAURICIO G. C. RESENDE

		Bias:
		Select elite solution and other solution.
		*/

        const individual *ind1, *ind2;
        if (pop.size() < 1)
        {
            cout << "Error - Pop size < 1";
            exit(1);
        }

        //select ind1 as random elite and ind2 as general random
        ind1 = ind2 = pop.select_bests(params.n_elite);

        while (ind2 == ind1)
        {
            if (params.second_parent_non_elite)
            {
                ind2 = pop.select_worst(pop.size() - params.n_elite);
            }
            else
            {
                ind2 = pop.select_worst(pop.size());
            }
        }

        return make_pair(ind1, ind2);
    }

    void apply_path_move(int &curr_dist, vvi &curr_indv, const vvi &target_indv, vector<pair<int, int>> &wrong_pos, vvi &f)
    {
        //f[mi][ji] -> position of ji in curr_indv
        bool move_applyed = false;

        while (!move_applyed && wrong_pos.size() > 0)
        {
            //choose random move
            int move_i = rand_gen() % wrong_pos.size();
            auto wpos = wrong_pos[move_i];
            //remove wrong position from moves
            swap(wrong_pos[move_i], wrong_pos[wrong_pos.size() - 1]);
            wrong_pos.pop_back();
            //========

            int mi = wpos.first;
            int pi = wpos.second;
            int ji = target_indv[mi][pi];
            int p2 = f[mi][ji];
            //job ji has to be in position pi of machine mi+1 (currently in p2)
            if (params.log_mode)
            {
                cout << "Job " << ji << " has to be in position " << pi
                     << " of machine " << (mi + 1) << " (currently in " << p2 << ")" << endl;
            }

            //one move can fix two positions, so need to test if its realy wrong pos.
            //lazy strategy
            if (p2 == pi)
                continue;

            move_applyed = true;
            curr_dist--;
            int j2 = curr_indv[mi][pi];
            if (params.log_mode)
            {
                //(ji from p2 to pi and j2 from pi to p2)
                cout << "Applying swap (" << j2 << "," << ji << ")" << endl;
            }

            if (target_indv[mi][p2] == j2)
            {
                //two fix move
                curr_dist--;
                //lazy strategy = not remove wrong pos now. ignore later.

                if (params.log_mode)
                {
                    cout << "Double fix move applyed." << endl;
                }
            }

            swap(curr_indv[mi][p2], curr_indv[mi][pi]);
            f[mi][ji] = pi;
            f[mi][j2] = p2;

            if (params.log_mode)
            {
                cout << "After swap move" << endl
                     << inst->mch_perm_to_str(curr_indv)
                     << endl;
            }
        }
    }

    int distance(const individual &p1, const individual &p2)
    {
        int mi, pi, ji;
        const vvi &mch_perm1 = p1.mch_perm;
        const vvi &mch_perm2 = p2.mch_perm;
        int n = inst->n, m = inst->m;
        int nop = inst->nop;

        int dist = 0;
        for (mi = 0; mi < m; mi++)
        {
            for (pi = 0; pi < n; pi++)
            {
                if (mch_perm1[mi][pi] != mch_perm2[mi][pi])
                {
                    dist++;
                }
            }
        }
        return dist;
    }

    void generate_childrens_path_relink(const individual &p1, const individual &p2, population &off)
    {
        /*
		Path Relinking strategy based on
		A Tabu Search/Path Relinking Algorithm to Solve the Job Shop Scheduling Problem
		Bo Peng, Zhipeng Lu, T.C.E. Cheng
		*/

        const vvi &mch_perm1 = p1.mch_perm;
        const vvi &mch_perm2 = p2.mch_perm;
        int n = inst->n, m = inst->m;
        int nop = inst->nop;

        //f[mi][ji] -> position in indv 1
        vvi &f = aux1;
        vector<pair<int, int>> wrong_pos;

        int mi, pi, ji;
        for (mi = 0; mi < m; mi++)
        {
            for (pi = 0; pi < n; pi++)
            {
                ji = mch_perm1[mi][pi];
                f[mi][ji] = pi;

                if (mch_perm1[mi][pi] != mch_perm2[mi][pi])
                {
                    wrong_pos.push_back(make_pair(mi, pi));
                }
            }
        }

        int dist = wrong_pos.size();
        stats.parents_dist.add(dist);

        vvi curr_indv = mch_perm1;
        int curr_dist = dist;
        int jump_size = max(dist / params.break_path_blocks, 2.0);
        int border_cut = max(dist / params.denominator_cut_pathrelink_borders, 1.0);
        int n_infeasible_local = 0, n_feasible_local = 0;

        solution s_best(*inst), s_eval(*inst);
        double best_eval = -1;

        params.ls.config(params.ls_params_path_eval);
        bool stop_path = false;

        if (params.log_mode)
        {
            cout << "Generating path relink childrens of :" << endl;
            cout << inst->mch_perm_to_str(mch_perm1);
            cout << endl
                 << " to " << endl;
            cout << inst->mch_perm_to_str(mch_perm2);
            cout << endl;
            cout << "Path distance = " << dist << endl;
            cout << "Jump size = " << jump_size << endl;
        }

        while (curr_dist > jump_size && !stop_path)
        {
            //advance jump_size path moves
            int old_dist = curr_dist;
            while (curr_dist > (old_dist - jump_size))
            {
                apply_path_move(curr_dist, curr_indv, mch_perm2, wrong_pos, f);
            }
            if (params.log_mode)
            {
                cout << "Stop Current Distance = " << curr_dist << endl;
            }
            if (params.cut_pathrelink_borders &&
                    (curr_dist > dist - border_cut) ||
                (curr_dist < border_cut))
            {
                //Apply border cut
                continue;
            }

            //evaluate curr solution
            bool is_feasible = s_eval.config(curr_indv);
            if (!is_feasible)
            {
                if (params.log_mode)
                {
                    cout << "Infeasible. Maybe Repair. " << endl;
                }
                stats.n_infeasible++;
                n_infeasible_local++;
                if (params.repair_infeasible)
                {
                    is_feasible = s_eval.repair();
                }
            }

            if (!is_feasible)
            {
                if (params.log_mode)
                {
                    cout << "Infeasible. No more repair." << endl;
                }
                continue;
            }
            n_feasible_local++;

            //Light TS to evaluate solution
            if (params.do_path_eval_ls)
            {
                params.ls.run_ls(s_eval);
            }

            double curr_eval = s_eval.fo_online;
            if (params.log_mode)
            {
                cout << "Feasible. Eval = " << curr_eval << endl;
            }
            //todo add something to take distance account

            //test if is the new best
            if (best_eval < 0 || curr_eval < best_eval)
            {
                s_best = s_eval;
                best_eval = curr_eval;
                if (params.log_mode)
                {
                    cout << "New best of path relink" << endl;
                }
            }
        }

        if (params.log_mode)
        {
            cout << "End of path crossover" << endl;
            cout << "#Infeaslible Sols = " << n_infeasible_local << endl;
            cout << "#Feasible Sols = " << n_feasible_local << endl;
        }

        //add best in offspring
        if (best_eval >= 0)
        {
            if (params.do_offspring_intensive_ls)
            {
                params.ls.config(params.ls_params_offspring);
                params.ls.run_ls(s_best);
            }

            if (params.log_mode)
            {
                cout << "Offspring generated. Makespan = " << s_best.fo_online << endl;
                cout << inst->mch_perm_to_str(s_best.my_mch_perm) << endl;
            }

            off.add(individual(s_best));
        }
    }

    void generate_childrens(const individual &elite, const individual &normal, population &off)
    {
        generate_childrens_path_relink(normal, elite, off);
        if (params.generate_both_directions)
        {
            generate_childrens_path_relink(elite, normal, off);
        }
    }

    void generate_offspring(population &pop, population &off)
    {
        off.clear();
        int i, n_generated = 0, to_generate = params.n_off_to_gen;

        int max_par_no_off = params.max_parents_no_offspring, n_par_no_off = 0;
        bool stop_try_gen_offspring = false;

        while (off.size() < to_generate && !stop_try_gen_offspring)
        {
            //select parents
            auto parents = select_parents(pop);
            //generate and add childrens
            auto old_size = off.size();
            generate_childrens(*parents.first, *parents.second, off);

            if (off.size() == old_size)
            {
                n_par_no_off++;
            }

            if (n_par_no_off >= max_par_no_off)
            {
                stop_try_gen_offspring = true;
            }
        }
        if (params.log_mode)
        {
            cout << "Total " << off.size() << "offspring generated." << endl;
            cout << "Mean parents distance = " << stats.parents_dist.get_mean() << endl;
        }
    }

    void update_pop(population &pop, population &off)
    {
        /*
		Methos based on
		A BIASED RANDOM-KEY GENETIC ALGORITHM FOR JOB-SHOP SCHEDULING
		JOS� FERNANDO GON�ALVES AND MAURICIO G. C. RESENDE

		The best n_elite are moved to next generation.
		The offspring is all included in new generation.
		The rest of new generation is completed with random solutions generated (mutation).
		*/

        int n_elite = params.n_elite;
        int pop_size_target = params.initial_pop_size;

        pop.keep_best(n_elite);
        pop.consume(off); //off is emptied

        while (pop.size() < pop_size_target)
        {
            //add random generated solution
            individual indv = gen_new_individual();
            pop.add(indv);
        }

        //keep pop sorted.
        pop.sort();
    }

    void run_mh()
    {
        /*
		Hybrid Genetic Algorithm
		1 - Generate Initial Population
		2 - Generate Offspring
		3 - Generate Mutations (Introducing Random Generated Solutions)
		4 - Update Population
		*/
        stats.start_time = get_time::now();
        // ========== MH Start ==============

        population pop, off;
        gen_init_pop(pop);
        update_best_fitness(pop);

        for (stats.curr_iteration = 0; !stop_mh(pop); stats.curr_iteration++)
        {
            auto start_it_time = get_time::now();
            //cout << "it " << stats.curr_iteration << endl;

            //===== MH ITERATION ====
            if (params.log_mode)
            {
                cout << " Begining of MH Iterations." << endl
                     << " Best = " << stats.best_fo_found << endl
                     << " Pop Size = " << pop.size()
                     << " #no improve its = " << stats.n_iter_no_global_improve << endl;
            }

            generate_offspring(pop, off);
            update_pop(pop, off);
            update_best_fitness(pop);
            //cout << "mean_distance = " << stats.parents_dist.get_mean() << endl;

            if (params.log_mode)
            {
                cout << " End of MH Iterations." << endl;
            }

            //========================
            auto end_it_time = get_time::now();
            auto it_time = to_time(end_it_time - start_it_time);
            stats.it_time_stream.add(it_time);
        }

        //Save the best solution found by the MH
        if (pop.size() == 0)
        {
            stats.best_sol_found = solution(*inst);
        }
        else
        {
            auto best_indv = pop.get_best();
            stats.best_sol_found.config(best_indv.mch_perm);
        }

        // ========== MH End ================

        auto end_mh_time = get_time::now();
        stats.total_time = to_time(end_mh_time - stats.start_time);
    }
};

// ================= Instance Selection =============================

class instance_selector
{
    /* Instance format:
	<header line>
	<n jobs> <n machines> ....
	Times
	<Time matrix: T_i_j = process time of operation j of job i.
	Machiens
	<Machine matrix: M_i_j = machine that operation j of job i uses.
	*/

    bool is_unique_inst = false, generated = false;
    string unique_inst_name;

public:
    string insts_folder_path = "./instances/";

    int ib, ie, i_cur;

    string inst_name, path, inst_prefix;

    instance_selector(
        int _ib = 1,
        int _ie = 10,
        string prefix = "Ta")
    {
        ib = _ib;
        ie = _ie;
        i_cur = ib - 1;
        inst_prefix = prefix;
    }

    instance_selector(string _unique_inst_name)
    {
        is_unique_inst = true;
        this->unique_inst_name = _unique_inst_name;
    }

private:
    string generate_path()
    {
        stringstream s;
        s << inst_prefix << setw(2) << setfill('0') << i_cur << ".txt";
        return insts_folder_path + s.str();
    }

public:
    string getNext()
    {
        if (is_unique_inst)
        {
            if (generated)
                return "";
            generated = true;

            return insts_folder_path + unique_inst_name;
        }
        i_cur++;

        if (i_cur > ie || i_cur < ib)
            return "";
        auto ret = generate_path();

        return ret;
    }

    string signature(char sep = ' ')
    {
        stringstream s;
        if (is_unique_inst)
        {
            s << unique_inst_name.substr(0, unique_inst_name.find('.'));
        }
        else
        {
            s << inst_prefix << setw(2) << setfill('0') << i_cur;
        }

        return s.str();
    }
};

// ========================= Experiments =======================

void do_test_intances_check()
{
    // ======= Instance selection config ===========
    int min_inst_idx = 1;
    int max_inst_idx = 3;
    // ======================

    instance_selector inst_selector(min_inst_idx, max_inst_idx);
    //instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);

        string solpath = inst_selector.insts_folder_path + inst_selector.signature() + ".sol";

        ifstream solfile;
        solfile.open(solpath);
        if (solfile.bad())
        {
            cout << "Bad file" << endl;
            return;
        }

        int n = inst.n;
        int m = inst.m;
        //build sol
        vvi mch_perm(m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int ji;
                solfile >> ji;
                mch_perm[i].push_back(ji + 1);
            }
        }

        solfile.close();

        cout << "Testing instance check." << endl;
        if (!inst.hard_check(mch_perm))
            cout << "Error Instance Check" << endl;

        int fo;
        inst.compute_fo(mch_perm, fo);
        cout << "Inst Makespan = " << fo << endl;

        cout << "Testing solution loading" << endl;
        solution s(inst);
        auto ok = s.config(mch_perm);
        if (!ok)
        {
            cout << "Error Solution Load" << endl;
            return;
        }
        cout << "Sol Makespan = " << s.fo_online << endl;

        if (s.fo_online != fo)
        {
            cout << "Test fail: objective functions dont match." << endl;
            return;
        }
    }
    cout << "End test" << endl;
}

void do_test_job_perm_loading()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);

        string solpath = inst_selector.insts_folder_path + inst_selector.signature() + ".sol2";
        ifstream solfile;
        solfile.open(solpath);
        if (solfile.bad())
        {
            cout << "Bad file" << endl;
            return;
        }

        int n = inst.n;
        int m = inst.m;
        //build sol - operation order
        vi job_perm(n * m);
        for (int i = 0; i < n * m; i++)
        {
            solfile >> job_perm[i];
        }
        solfile.close();
        //job perm => op perm
        inst.job_perm_to_op_perm(job_perm);

        cout << "Testing solutiono order loading" << endl;

        cout << "Loading operation order: ";
        for (int i = 0; i < job_perm.size(); i++)
            cout << pii_to_str(inst.pair_job_pos(job_perm[i])) << " ";
        cout << endl;

        solution s(inst);
        bool ok = s.config(job_perm);

        if (!ok)
        {
            cout << "Error Solution Load" << endl;
        }
        cout << "Sol Makespan = " << s.fo_online << endl;
    }
    cout << "End test" << endl;
}

void do_test_sol_gen()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        inst.print();

        cout << "Testing solution generation" << endl;
        solution_generator sol_gen(inst);
        int nit = 100;
        for (int i = 0; i < nit; i++)
        {
            solution s = sol_gen.generate();
            cout << "Solution " << i << " generated:" << endl;
            s.print(true);
        }
    }
    cout << "End test" << endl;
}

void do_simple_test_ls()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        inst.print();

        cout << "Testing basic ls" << endl;

        string solpath = inst_selector.insts_folder_path + inst_selector.signature() + ".sol2";
        ifstream solfile;
        solfile.open(solpath);
        if (solfile.bad())
        {
            cout << "Bad file" << endl;
            return;
        }

        int n = inst.n;
        int m = inst.m;
        //build sol - operation order
        vi job_perm(n * m);
        for (int i = 0; i < n * m; i++)
        {
            solfile >> job_perm[i];
        }
        solfile.close();
        //job perm => op perm
        inst.job_perm_to_op_perm(job_perm);

        cout << "Loading operation order: ";
        for (int i = 0; i < job_perm.size(); i++)
            cout << pii_to_str(inst.pair_job_pos(job_perm[i])) << " ";
        cout << endl;

        solution s(inst);
        s.config(job_perm);

        s.print();

        ls_engine ls;
        ls.config();
        ls.params.debug_mode = true;

        cout << "Starting LS" << endl;
        //restart_global_clock(2, true); //2 sec
        ls.run_ls(s);
        cout << "End of LS" << endl;

        cout << "Hard check test" << endl;
        if (!inst.hard_check(ls.best_sol_found(), ls.best_fo_found()))
        {
            cout << "Failed hard check." << endl;
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Sol. Found = " << ls.best_fo_found() << endl;
        }
    }
    cout << "End test" << endl;
}

void do_hard_correctness_test_ls()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector(1, 40);

    string inst_path = "";
    cout << "Testing Heavy Ls" << endl;

    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        solution_generator sol_gen(inst);

        cout << "Instance " << inst_selector.signature() << endl;
        int n_try = 100;

        ls_engine ls;
        ls.config();
        ls.params.max_n_moves = 20000;
        ls.params.max_moves_no_global_improve = 1000;

        solution s(inst);
        int min_makespan = I32_INF;
        stats_stream<long> tstream("time[s]");

        for (int it = 1; it <= n_try; it++)
        {
            //cout << "It = " << it << endl;
            s = sol_gen.generate();
            ls.run_ls(s);

            //cout << "hard check test" << endl;
            if (!inst.hard_check(ls.best_sol_found(), ls.best_fo_found()))
            {
                cout << "Failed hard check." << endl;
                exit(1);
            }
            else
            {
                //cout << "Hard check ok. Makespan = " << s->fo_online << endl;
                min_makespan = min(min_makespan, ls.best_fo_found());
                tstream.add(ls.stats.total_time);
            }
        }
        cout << "Min Makespaen Found = " << min_makespan
             << " (" << tstream.get_mean() << " [s])"
             << endl;
    }
    cout << "End test" << endl;
}

void do_time_test_ls()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector(1, 20);

    string inst_path = "";
    cout << "Testing Heavy Ls" << endl;

    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        //inst.print();
        solution_generator sol_gen(inst);

        cout << "Instance " << inst_selector.signature() << endl;
        int n_try = 100;

        ls_engine ls;
        ls.config();
        ls.params.max_n_moves = 10000;
        ls.params.max_moves_no_global_improve = 1000;
        ls.params.tabu_list_factor = 2;

        solution s(inst);
        int min_makespan = I32_INF;
        stats_stream<long> tstream("time[s]");

        for (int it = 1; it <= n_try; it++)
        {
            //cout << "It = " << it << endl;
            s = sol_gen.generate();
            ls.run_ls(s);

            min_makespan = min(min_makespan, ls.best_fo_found());
            tstream.add(ls.stats.total_time / 1000);
        }

        cout << "Min Makespaen Found = " << min_makespan
             << " (" << tstream.get_mean() << " [s])"
             << endl;
    }
    cout << "End test" << endl;
}

string default_header = "Instance BestSol AvgSol AvgTime NIt";
void do_inst_ls_solve()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector(1, 40);

    string inst_path = "";
    cout << default_header << endl;

    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        solution_generator sol_gen(inst);

        cout << inst_selector.signature();
        int n_try = 100;

        ls_engine ls;
        ls.config();

        solution s(inst);
        long min_makespan = I32_INF, n_stop_no_move = 0;
        stats_stream<double> time_stream("time[s]");
        stats_stream<long> fo_stream("Makespan");

        for (int it = 1; it <= n_try; it++)
        {
            s = sol_gen.generate();
            ls.run_ls(s);
            if (ls.stats.stop_because_no_move)
                n_stop_no_move++;

            fo_stream.add(ls.best_fo_found());
            time_stream.add(ls.stats.total_time / 1000.0);
            min_makespan = fo_stream.get_min();
            if (min_makespan == ls.best_fo_found())
            {
                if (!inst.hard_check(ls.best_sol_found(), min_makespan))
                {
                    cout << "Error" << endl;
                    exit(1);
                }
            }
        }

        cout << " " << min_makespan
             << " " << fo_stream.get_mean()
             << " " << time_stream.get_mean()
             << " " << (n_stop_no_move)
             << endl;
    }
}

void do_simple_test_mh()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        inst.print();

        cout << "Testing (simple) MH" << endl;

        mh_engine mh;
        mh.config();
        mh.params.initial_pop_size = 4;
        mh.params.curr_pop_size = 4;
        mh.params.n_elite = 1;
        mh.params.n_off_to_gen = 2;

        mh.params.log_mode = true;

        solution s = mh.run(inst);
        if (s.fo_online != mh.get_best_fo())
        {
            cout << "Error s.fo_online != mh.get_best_fo() ";
            exit(1);
        }

        cout << "Hard check test" << endl;
        if (!inst.hard_check(s.my_mch_perm, s.fo_online))
        {
            cout << "Failed hard check." << endl;
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Makespan = " << mh.get_best_fo() << endl;
        }
    }
    cout << "End test" << endl;
}

void do_hard_test_mh()
{
    // ======= Instance selection config ===========

    //instance_selector inst_selector(min_inst_idx, max_inst_idx);
    instance_selector inst_selector(1, 10);

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        inst.print();

        cout << "Testing (hard) MH" << endl;

        mh_engine mh;
        mh.config();

        solution s = mh.run(inst);
        if (s.fo_online != mh.get_best_fo())
        {
            cout << "Error s.fo_online != mh.get_best_fo() ";
            exit(1);
        }

        cout << "Hard check test" << endl;
        if (!inst.hard_check(s.my_mch_perm, s.fo_online))
        {
            cout << "Failed hard check." << endl;
            exit(1);
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Makespan = " << mh.get_best_fo() << endl;
        }
    }
    cout << "End test" << endl;
}

void do_inst_mh_solve()
{

    // instance selection config
    instance_selector inst_selector(1, 40); //Ta01 - Ta40

    string inst_path = "";
    cout << default_header << endl;
    int n_try = 10;
    cout << "N_Runs = " << n_try << endl;

    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;

        instance inst;
        inst.initialize(inst_path);
        mh_engine mh;
        mh.config();

        cout << inst_selector.signature();

        long min_makespan = I32_INF;
        stats_stream<double> time_stream("time[s]");
        stats_stream<long> fo_stream("Makespan");
        stats_stream<int> n_it_stream("N IT");

        for (int it = 1; it <= n_try; it++)
        {
            solution s = mh.run(inst);
            auto old_makespan = min_makespan;

            fo_stream.add(mh.get_best_fo());
            time_stream.add(mh.stats.total_time / 1000.0);
            min_makespan = fo_stream.get_min();
            n_it_stream.add(mh.stats.curr_iteration);

            if (min_makespan < old_makespan)
            {
                if (!inst.hard_check(s.my_mch_perm, min_makespan))
                {
                    cout << "Error" << endl;
                    exit(1);
                }
            }
        }

        cout << " " << min_makespan
             << " " << fo_stream.get_mean()
             << " " << time_stream.get_mean()
             << " " << n_it_stream.get_mean()
             << endl;
    }
}

// ==============================================================

int main()
{
    
    // =========== IO config =================
    string output_path = "output/jssp_out.txt";
    std::ofstream out(output_path);

    if (out)
    {
        cout << "Output file: " << output_path << endl;
    }
    else
    {
        cout << "Unable to create output file: " << output_path << endl;
        exit(1);
    }

    auto coutbuf = std::cout.rdbuf(out.rdbuf());

    //do_test_intances_check();
    //do_test_job_perm_loading();
    //do_test_sol_gen();
    //do_simple_test_ls();
    //do_hard_correctness_test_ls();
    //do_time_test_ls();
    //do_inst_ls_solve();
    do_simple_test_mh();
    //do_hard_test_mh();
    //do_inst_mh_solve();

    std::cout.rdbuf(coutbuf);
}
