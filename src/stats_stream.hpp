#pragma once
#include "global.hpp"
#include <list>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>

template <typename Tnum>
class stats_stream
{
    using string = std::string;
    
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

    std::pair<Tnum, Tnum> get_min_max();

    double get_sd();

    Tnum get_min();
    
    Tnum get_max();

    bool empty();

    void print(bool print_total = false);

};

template <typename T>
stats_stream<T>::stats_stream() {}

template <typename T>
stats_stream<T>::stats_stream(std::string _name) : name(_name) {}

template <typename T>
stats_stream<T>::stats_stream(std::string _name, bool _do_sd, bool _save_seq, bool _over_time)
    : name(_name), do_sd(_do_sd), save_seq(_save_seq), is_over_time(_over_time){};

template <typename T>
void stats_stream<T>::clear()
{
    firstSample = true;
    total = 0, total_2 = 0, nsample = 0;
    seq.clear();
    time_seq.clear();
}

template <typename T>
void stats_stream<T>::add(T x, double time_point)
{
    nsample++;
    firstSample ? (minv = x) : (minv = std::min(minv, x));
    firstSample ? (maxv = x) : (maxv = std::max(maxv, x));
    firstSample = false;
    total += x;
    if (total > overflow_test)
        std::cout << "Error Overflow risk" << std::endl, exit(1);
    if (do_sd)
    {
        total_2 += x * x;
        if (total_2 > overflow_test)
        {
            std::cout << "Error Overflow risk" << std::endl;
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


template <typename T>
T stats_stream<T>::get_total() { return total; }

template <typename T>
double stats_stream<T>::get_mean() { return nsample != 0 ? ((double)total / nsample) : 0; }

template <typename T>
std::pair<T, T> stats_stream<T>::get_min_max() { return make_pair(minv, maxv); }

template <typename T>
double stats_stream<T>::get_sd()
{
    if (nsample == 0)
        return 0;
    //sqrt(E[X^2] - E[X]^2)
    auto var = (total_2 / nsample) - ((long double)total / nsample) * ((double)total / nsample);
    return sqrtl(var);
}

template <typename T>
T stats_stream<T>::get_min() { return minv; };

template <typename T>
T stats_stream<T>::get_max() { return maxv; };

template <typename T>
bool stats_stream<T>::empty() { return nsample == 0; }

template <typename T>
void stats_stream<T>::print(bool print_total)
{
    std::cout << name << ": "
            << "min = " << minv << " max = " << maxv << " mean = " << get_mean();
    if (do_sd)
        std::cout << " sd = " << get_sd();
    if (print_total)
        std::cout << " total = " << total;
    std::cout << std::endl;
}