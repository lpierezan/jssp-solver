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
#include "operation.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include "solution_generator.hpp"
#include "ls_engine.hpp"

using namespace std;

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
                ret << p->to_string(inst) << endl;
            return ret.str();
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
        int mi, pi;
        const vvi &mch_perm1 = p1.mch_perm;
        const vvi &mch_perm2 = p2.mch_perm;
        int n = inst->n, m = inst->m;
        
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
                    ((curr_dist > dist - border_cut) || (curr_dist < border_cut))
                    )
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
        int to_generate = params.n_off_to_gen;

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

    string signature()
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
