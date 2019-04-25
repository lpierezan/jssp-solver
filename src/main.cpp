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
#include "mh_engine.hpp"
#include "instance_selector.hpp"


using namespace std;

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
