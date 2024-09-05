#include "Utilities.hpp"

namespace PSF
{

int factorial(int n) {
  int output = 1;
  for(int i = 1; i <= n; i++) {
    output *= i;
  }
  return output;
}

// Recursive function to find partitions of an integer
void findPartitions(int n, std::vector<int>& current, std::vector<std::vector<int>>& result, int start) {
  if (n == 0) {
    // When we have successfully partitioned the integer, add to result
    result.push_back(current);
    return;
  }

  // Try numbers from `start` to `n`
  for (int i = start; i <= n; ++i) {
    // Add this part to the current partition
    current.push_back(i);

    // Recursively partition the remaining sum `n - i`
    findPartitions(n - i, current, result, i);

    // Backtrack to remove the last added element
    current.pop_back();
  }
}

// Function to get all partitions of a given integer
std::vector<std::vector<int>> getPartitions(int n) {
  std::vector<std::vector<int>> result;
  std::vector<int> current;

  // Start partitioning with the smallest part being 1
  findPartitions(n, current, result, 1);

  return result;
}

// Recursive function to generate permutations
void generatePermutations(std::vector<int>& arr, int start, std::vector<std::vector<int>>& result) {
  // Base case: If we've reached the end of the array, add the permutation to the result
  if (start == arr.size()) {
    result.push_back(arr);
    return;
  }

  // Iterate over the elements and swap them to generate permutations
  for (int i = start; i < arr.size(); ++i) {
    std::swap(arr[start], arr[i]);
    generatePermutations(arr, start + 1, result); // Recurse with the next start
    std::swap(arr[start], arr[i]); // Backtrack to restore original state
  }
}

// Function to get all permutations of a given vector
std::vector<std::vector<int>> getPermutations(const std::vector<int>& vec) {
  std::vector<std::vector<int>> result;
  std::vector<int> arr = vec;

  // Start permutation generation from the first index
  generatePermutations(arr, 0, result);
  removeDuplicates(result);
  return result;
}


// A custom comparator for a set of vectors
struct VectorCompare {
  bool operator()(const std::vector<int>& a, const std::vector<int>& b) const {
    // Compares two vectors by lexicographic order
    return a < b;
  }
};

// Function to remove duplicates from a vector of vectors
void removeDuplicates(std::vector<std::vector<int>>& vecOfVecs) {
  // Create a set to store unique vectors using the custom comparator
  std::set<std::vector<int>, VectorCompare> uniqueVectors;

  // Insert all vectors into the set to remove duplicates
  for (const auto& vec : vecOfVecs) {
    uniqueVectors.insert(vec);
  }

  // Clear the original vector and copy the unique vectors back to it
  vecOfVecs.clear();
  vecOfVecs.assign(uniqueVectors.begin(), uniqueVectors.end());
}

// Function to recursively generate subsets
void generateSubsetsHelper(const std::vector<int>& original, std::vector<int>& currentSubset, std::vector<std::vector<int>>& subsets, int start, int n) {
  // If the current subset has reached the desired length, add it to the subsets list
  if (currentSubset.size() == n) {
    subsets.push_back(currentSubset);
    return;
  }

  // Iterate through the remaining elements and recursively generate subsets
  for (int i = start; i < original.size(); ++i) {
    // Add the current element to the subset
    currentSubset.push_back(original[i]);

    // Recursively generate the rest of the subset
    generateSubsetsHelper(original, currentSubset, subsets, i + 1, n);

    // Backtrack: remove the last element to explore other subsets
    currentSubset.pop_back();
  }
}

// Wrapper function to generate subsets of length n
std::vector<std::vector<int>> generateSubsets(const std::vector<int>& original, int n) {
  std::vector<std::vector<int>> subsets;
  std::vector<int> currentSubset;

  if (n > original.size() || n < 0) {
    return subsets; // Invalid case
  }

  // Sort original vector (optional, but often desirable)
  std::vector<int> sortedOriginal = original;
  sort(sortedOriginal.begin(), sortedOriginal.end());

  // Start the recursive process
  generateSubsetsHelper(sortedOriginal, currentSubset, subsets, 0, n);

  return subsets;
}

}