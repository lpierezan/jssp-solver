#pragma once

#include <random>
#include "global.hpp"

class permutation
{
public:
    static bool check(vi v, int min = 1, int max = -1);

    static vi generate_rand(int min, int max, std::default_random_engine &rand_gen);
};