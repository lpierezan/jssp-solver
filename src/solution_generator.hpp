#pragma once
#include "instance.hpp"
#include "solution.hpp"
#include <random>

class solution_generator
{
    instance *my_inst;
    std::default_random_engine my_rand_gen = std::default_random_engine(0);

public:
    solution_generator();
    
    solution_generator(instance &inst, std::default_random_engine _rand_gen = std::default_random_engine(0));

    void reset_rand(int seed = 0);

    solution generate();
    
};