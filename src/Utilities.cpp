#include "Utilities.hpp"

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

// Recursive function to generate permutations
void generatePermutations(std::vector<int>& arr, int start, const std::vector<bool>& flavor, std::vector<std::vector<int>>& result) {
  // Base case: If we've reached the end of the array, add the permutation to the result
  if (start == arr.size()) {
    result.push_back(arr);
    return;
  }

  // Iterate over the elements and swap them to generate permutations
  for (int i = start; i < arr.size(); ++i) {
    if(flavor[start] != flavor[i]) {
      std::swap(arr[start], arr[i]);
      generatePermutations(arr, start + 1, flavor, result); // Recurse with the next start
      std::swap(arr[start], arr[i]); // Backtrack to restore original state
    }
    else{
      result.push_back(arr);
    }
  }
}

// Function to get all permutations of a given vector
std::vector<std::vector<int>> getPermutations(const std::vector<int>& vec, const std::vector<bool>& flavor) {
  std::vector<std::vector<int>> result;
  std::vector<int> arr = vec;

  // Start permutation generation from the first index
  generatePermutations(arr, 0, flavor, result);
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
