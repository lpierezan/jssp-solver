/*
TODO Improve to some test framework.
*/

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

#include <iostream>
#include <fstream>

using namespace std;

bool do_test_intances_check()
{
    // ======= Instance selection config ===========
    int min_inst_idx = 1;
    int max_inst_idx = 5;
    // ======================

    instance_selector inst_selector(min_inst_idx, max_inst_idx);
    //instance_selector inst_selector("test.txt");

    string inst_path = "";
    while (true)
    {
        inst_path = inst_selector.getNext();
        if (inst_path == "")
            break;
        
        cout << "Testing instance: " << inst_path << endl;

        instance inst;
        inst.initialize(inst_path);

        string solpath = inst_selector.insts_folder_path + inst_selector.signature() + ".sol";

        ifstream solfile;
        solfile.open(solpath);
        if (solfile.bad())
        {
            cout << "Bad solution file" << endl;
            return false;
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

        cout << "Instance hard check...." << endl;
        if (!inst.hard_check(mch_perm))
        {
            cout << "Hard Check Failed." << endl;
            return false;
        }
            

        cout << "Compute objective function..." << endl;
        int fo;
        inst.compute_fo(mch_perm, fo);
        cout << "Instance Makespan = " << fo << endl;

        cout << "Solution loading..." << endl;
        solution s(inst);
        auto ok = s.config(mch_perm);
        if (!ok)
        {
            cout << "Error Solution Load" << endl;
            return false;
        }
        cout << "Solution Makespan = " << s.fo_online << endl;

        cout << "Compare objective functions..." << endl;

        if (s.fo_online != fo)
        {
            cout << "Test fail: objective functions dont match." << endl;
            return false;
        }
    }

    cout << "end test" << endl;

    return true;
}

bool do_test_job_perm_loading()
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
            cout << "Solution bad file." << endl;
            return false;
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

        cout << "Testing solution operations order loading..." << endl;

        cout << "Loading operation order: ";
        for (int i = 0; i < job_perm.size(); i++)
            cout << pii_to_str(inst.pair_job_pos(job_perm[i])) << " ";
        cout << endl;

        solution s(inst);
        bool ok = s.config(job_perm);

        if (!ok)
        {
            cout << "Error Solution Load" << endl;
            return false;
        }

        cout << "Sol Makespan = " << s.fo_online << endl;
    }

    cout << "end test" << endl;

    return true;
}

bool do_test_sol_gen()
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
    return true;
}

bool do_simple_test_ls()
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
            cout << "Solution bad file" << endl;
            return false;
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
            return false;
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Sol. Found = " << ls.best_fo_found() << endl;
        }
    }

    cout << "End test" << endl;
    return true;
}

bool do_hard_correctness_test_ls()
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
                return false;
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
    return true;
}

bool do_time_test_ls()
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
    return true;
}

bool do_simple_test_mh()
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
            return false;
        }

        cout << "Hard check test" << endl;
        if (!inst.hard_check(s.my_mch_perm, s.fo_online))
        {
            cout << "Failed hard check." << endl;
            return false;
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Makespan = " << mh.get_best_fo() << endl;
        }
    }

    cout << "End test" << endl;
    return true;
}

bool do_hard_test_mh()
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
            return false;
        }

        cout << "Hard check test" << endl;
        if (!inst.hard_check(s.my_mch_perm, s.fo_online))
        {
            cout << "Failed hard check." << endl;
            return false;
        }
        else
        {
            cout << "Hard check ok." << endl;
            cout << "Best Makespan = " << mh.get_best_fo() << endl;
        }
    }
    cout << "End test" << endl;
    return true;
}

int main()
{

    cout << "Starting Instance tests..." << endl;
    int n_fail = 0;

    if(!do_test_intances_check()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_test_job_perm_loading()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_test_sol_gen()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_simple_test_ls()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_hard_correctness_test_ls()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_time_test_ls()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_simple_test_mh()){
        n_fail++;
        cout << "Failed test!" << endl;
    }

    if(!do_hard_test_mh()){
        n_fail++;
        cout << "Failed test!" << endl;
    }
    
    if(n_fail == 0)
    {
        cout << "Passed :)" << endl;
    }else{
        cout << "Failed " << n_fail << " tests. :(" << endl;
    }
    
    return 0;
}