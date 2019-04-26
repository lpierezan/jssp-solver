#include "global.hpp"
#include "experiments.hpp"
#include "instance.hpp"
#include "instance_selector.hpp"
#include "mh_engine.hpp"
#include "ls_engine.hpp"
#include "solution.hpp"

namespace experiments
{

using namespace std;

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

}