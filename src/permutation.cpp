#include <vector>
#include <random>
#include <algorithm>
#include "global.hpp"
#include "permutation.hpp"

bool permutation::check(vi v, int min, int max)
{
    //[min , max]
    if (max == -1)
        max = v.size();
    if (max - min + 1 != v.size())
        return false;
    std::vector<bool> mark(v.size(), false);
    for (int i = 0; i < v.size(); i++)
    {
        int p = v[i] - min;
        if (p < 0 || p >= mark.size())
            return false;
        if (mark[p])
            return false;

        mark[p] = true;
    }
    return true;
}

vi permutation::generate_rand(int min, int max, std::default_random_engine &rand_gen)
{
    vi perm(max - min + 1);
    for (int i = 0; i < perm.size(); i++)
    {
        perm[i] = i + min;
    }
    shuffle(ALL(perm), rand_gen);
    return perm;
}