#include "global.hpp"
#include "ls_engine.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include "stats_stream.hpp"
#include "time_control.hpp"
#include <vector>
#include <unordered_set>
#include <deque>

using namespace std;

string ls_engine::move::to_string(instance *inst)
{
    return pii_to_str(inst->pair_job_pos(op1)) + " " + (fw ? "fw" : "bw") + " " +
            pii_to_str(inst->pair_job_pos(op2));
}

bool ls_engine::move::check(solution &s)
{
    auto p1 = s.my_inst->pair_job_mch(op1);
    auto p2 = s.my_inst->pair_job_mch(op2);

    if (p1.first == p2.first)
        return false;
    if (p1.second != p2.second)
        return false;

    return true;
}

void ls_engine::move_mem_data::init(bool _fw)
{
    fw = _fw;
    int x_fw = fw ? 2 : 1;
    seed ^= hsh_int(x_fw) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

bool ls_engine::move_mem_data::eqp_cmp::operator()(const move_mem_data &move_mem1, const move_mem_data &move_mem2) const
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

void ls_engine::move_mem_data::add_value(int x)
{
    //compund hash similar to boost implementation.
    seed ^= hsh_int(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    block_idxs.push_back(x);
}

size_t ls_engine::move_mem_data::move_hash::operator()(const move_mem_data &move) const
{
    return move.seed;
}

string ls_engine::move_mem_data::to_string()
{
    string ret = (fw ? "fw " : "bw ");
    for (int x : block_idxs)
    {
        ret += std::to_string(x) + " ";
    }
    return ret;
}

void ls_engine::ls_stats::reset(solution &s)
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

void ls_engine::ls_parameters::compute(solution &s)
{
    tabu_list_size = s.n * s.m * (tabu_list_factor >= 0 ? tabu_list_factor : 1);
}

ls_engine::ls_parameters::ls_parameters() {}

ls_engine::ls_engine()
{
}

//Parameters Input
void ls_engine::config(ls_parameters _params)
{
    params = _params;
}

//String signature
string ls_engine::signature()
{
    //TODO
    string ret = "";
    return ret;
}

void ls_engine::run_ls(solution &initial_sol)
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

int ls_engine::best_fo_found()
{
    return stats.best_fo_found;
}

//returns the job permutation in each machine
vvi ls_engine::best_sol_found()
{
    //not neccessarly during the local search
    return stats.best_sol_found;
}

bool ls_engine::tabu_acc(move_mem_data &move_tabu_data, int eval)
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

void ls_engine::tabu_add(move_mem_data &move_tabu_data)
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

void ls_engine::create_tabu_move_data(move &move, vi &ms, vi &mp,
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

void ls_engine::tabu_list_size_update()
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

bool ls_engine::test_acyclic_move(int u, int v, bool is_fw)
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

void ls_engine::explore_block_moves(vi &block_idxs,
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

        int j;
        for (j = i - 1; j >= 1; j--)
        {
            int a = block_idxs[j];
            int asuc = js[a];

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

bool ls_engine::search_block(int block_begin, move &best_move, int &best_eval)
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

bool ls_engine::search_move(move &selected_move)
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

void ls_engine::update_tsv(int u, int v, bool fw,
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

void ls_engine::apply_move(move &move_apply)
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

bool ls_engine::stop_criterion_test()
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

void ls_engine::keep_best_sol_test()
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

void ls_engine::ls_start(solution &s)
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