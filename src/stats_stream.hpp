#pragma once
#include "global.hpp"
#include <list>

using std::string;
using std::pair;

template <typename Tnum>
class stats_stream
{
    /* Collect stream of data and compute online statistics (good when not storing sequence)
	*/
    string name;
    bool do_sd = false, save_seq = false,
         is_over_time = false;

    bool firstSample = true;
    Tnum total = 0, minv, maxv;
    long double total_2;

    long nsample = 0;
    std::list<Tnum> seq;
    std::list<double> time_seq;

    static const long long int overflow_test = ((long long)1 << 62);

public:
    stats_stream();
    
    stats_stream(string _name);
    
    stats_stream(string _name, bool _do_sd, bool _save_seq = false, bool _over_time = false);
    
    void clear();

    void add(Tnum x, double time_point = 0);

    Tnum get_total();

    double get_mean();

    pair<Tnum, Tnum> get_min_max();

    double get_sd();

    Tnum get_min();
    
    Tnum get_max();

    bool empty();

    void print(bool print_total = false);

};