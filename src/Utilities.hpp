#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

namespace PSF
{
/**
 * @brief Factorial function (n!)
 *
 * @param n Argument of the factorial function (n!)
 *
 */
int factorial(int n);

/**
 * @brief Generates all possible paritions of the integer n.
 *
 * Returns a list of all possible partitions.
 *
 * Example usage:
 * @code{.cpp}
 * getPartitions(4) = {{1,1,1,1},{1,1,2},{1,3},{2,2},{4}}
 * @endcode
 *
 * @param n Integer to generate partitions of.
 *
 */
std::vector<std::vector<int>> getPartitions(int n);

/**
 * @brief Generates all possible Permutations of a vector.
 *
 *
 * Example usage:
 * @code{.cpp}
 * getPermutations({1,2,3}) = {{1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1}}
 * @endcode
 *
 * @param vec Input vector to generate permutations of.
 *
 */
std::vector<std::vector<int>> getPermutations(const std::vector<int>& vec);

/**
 * @brief Removes duplicates of vectors from a list of vectors.
 *
 * @param vecOfVecs List of Vectors to remove duplicates from.
 *
 */
void removeDuplicates(std::vector<std::vector<int>>& vecOfVecs);

/**
 * @brief Generates all possible subsets (unordered) of a set with length n.
 *
 * @param original Original set to generate subsets from.
 * @param n Number of elements in the subset.
 *
 */
std::vector<std::vector<int>> generateSubsets(const std::vector<int>& original, int n);

}
#endif