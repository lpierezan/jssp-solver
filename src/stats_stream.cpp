#include "global.hpp"
#include "stats_stream.hpp"


using namespace std;

stats_stream::stats_stream() {}

stats_stream::stats_stream(string _name) : name(_name) {}

stats_stream::stats_stream(string _name, bool _do_sd, bool _save_seq = false, bool _over_time = false)
    : name(_name), do_sd(_do_sd), save_seq(_save_seq), is_over_time(_over_time){};

void stats_stream::clear()
{
    firstSample = true;
    total = 0, total_2 = 0, nsample = 0;
    seq.clear();
    time_seq.clear();
}

void stats_stream::add(Tnum x, double time_point = 0)
{
    nsample++;
    firstSample ? (minv = x) : (minv = min(minv, x));
    firstSample ? (maxv = x) : (maxv = max(maxv, x));
    firstSample = false;
    total += x;
    if (total > overflow_test)
        cout << "Error Overflow risk" << endl, exit(1);
    if (do_sd)
    {
        total_2 += x * x;
        if (total_2 > overflow_test)
        {
            cout << "Error Overflow risk" << endl;
            exit(1);
        }
    }
    if (save_seq)
    {
        seq.push_back(x);
        if (is_over_time)
            time_seq.push_back(time_point);
    }
}

Tnum stats_stream::get_total() { return total; }

double stats_stream::get_mean() { return nsample != 0 ? ((double)total / nsample) : 0; }

pair<Tnum, Tnum> stats_stream::get_min_max() { return make_pair(minv, maxv); }

double stats_stream::get_sd()
{
    if (nsample == 0)
        return 0;
    //sqrt(E[X^2] - E[X]^2)
    auto var = (total_2 / nsample) - ((long double)total / nsample) * ((double)total / nsample);
    return sqrtl(var);
}

Tnum stats_stream::get_min() { return minv; };
Tnum stats_stream::get_max() { return maxv; };
bool stats_stream::empty() { return nsample == 0; }

void stats_stream::print(bool print_total = false)
{
    cout << name << ": "
            << "min = " << minv << " max = " << maxv << " mean = " << get_mean();
    if (do_sd)
        cout << " sd = " << get_sd();
    if (print_total)
        cout << " total = " << total;
    cout << endl;
}