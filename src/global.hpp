#pragma once

#include <vector>
#include <random>
#include <chrono>

typedef std::vector<std::vector<int>> vvi;
typedef std::vector<int> vi;
#define ALL(v) (v).begin(), (v).end()
#define D_INF (DBL_MAX)
#define I32_INF (1 << 30)

extern std::default_random_engine rand_gen;

void reset_rand(int seed);

//clever online shuffle algorithm
int select_shuffle(vi &v, int &max_idx);

std::string pii_to_str(std::pair<int, int> p);
