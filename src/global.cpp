#include "global.hpp"

std::default_random_engine rand_gen = std::default_random_engine(0);

void reset_rand(int seed)
{
    rand_gen = std::default_random_engine(seed);
}

//clever online shuffle algorithm
int select_shuffle(vi &v, int &max_idx)
{
    auto pos = max_idx > 0 ? rand_gen() % (max_idx + 1) : 0;
    std::swap(v[pos], v[max_idx]);
    max_idx--;
    return v[max_idx + 1];
}

std::string pii_to_str(std::pair<int, int> p)
{ 
    return "(" + std::to_string(p.first) + "," + std::to_string(p.second) + ")";
}