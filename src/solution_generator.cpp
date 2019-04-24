#include "global.hpp"
#include "solution_generator.hpp"
#include <algorithm>
#include <iostream>

solution_generator::solution_generator() {}

solution_generator::solution_generator(instance &inst, std::default_random_engine _rand_gen)
{
    my_inst = &inst;
    my_rand_gen = _rand_gen;
}

void solution_generator::reset_rand(int seed)
{
    my_rand_gen = std::default_random_engine(seed);
}

solution solution_generator::generate()
{
    solution s(*my_inst);

    // == Permutation of all operations (Using repeted job_index permutation) =====
    vi job_rep_perm;
    for (int ji = 1; ji <= my_inst->n; ji++)
    {
        for (int pi = 1; pi <= my_inst->m; pi++)
        {
            job_rep_perm.push_back(ji);
        }
    }
    //Shuffle
    std::shuffle(ALL(job_rep_perm), my_rand_gen);

    //From job index => operation index
    my_inst->job_perm_to_op_perm(job_rep_perm);

    // ======= Config Solution ===========
    auto ok = s.config(job_rep_perm);
    if (!ok)
    {
        std::cout << "Error - Invalid solution generated. (not supposed to happen)" << std::endl;
        exit(0);
    }

    return s;
}