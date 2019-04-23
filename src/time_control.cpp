#include "global.hpp"
#include "time_control.hpp"

std::chrono::steady_clock::time_point global_time = get_time::now();
long long global_tl_sec = 0;
bool active_time_limit = false;

void restart_global_clock(long long tl_sec, bool active_tl_test)
{
    global_time = get_time::now();
    global_tl_sec = tl_sec;
    active_time_limit = active_tl_test;
}

bool tl_clock_check(long long &duration_ms)
{
    duration_ms = to_time(get_time::now() - global_time);
    return active_time_limit && ((duration_ms / 1000.0) >= global_tl_sec);
}
