#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

int factorial(int n);

std::vector<std::vector<int>> getPartitions(int n);

std::vector<std::vector<int>> getPermutations(const std::vector<int>& vec);

void removeDuplicates(std::vector<std::vector<int>>& vecOfVecs);


#endif