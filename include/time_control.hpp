#pragma once
#include "global.hpp"

//TODO MAKE A CLASS

// =============== Time/Clock Measurement ====================
using get_time = std::chrono::steady_clock;
#define to_time(x) (std::chrono::duration_cast<std::chrono::milliseconds>(x).count())

extern std::chrono::steady_clock::time_point global_time;
extern long long global_tl_sec;
extern bool active_time_limit;

void restart_global_clock(long long tl_sec, bool active_tl_test);

bool tl_clock_check(long long &duration_ms);