#ifndef TREE_HPP
#define TREE_HPP

#include <iostream>
#include <vector>
#include <string>

template <typename T>
struct TreeNode {
  T data;
  int level = 0;
  int index = 0;
  std::vector<TreeNode<T>*> children;
  TreeNode* parent;

  TreeNode(const T& value) : data(value) {}
};

template <typename T>
class Tree {
private:
  TreeNode<T>* root;
  void destroyTree(TreeNode<T>* node) {
    if (node != nullptr) {
      for (TreeNode<T>* child : node->children) {
        destroyTree(child);
      }
      delete node;
    }
  }

  // Recursive function to copy a tree
  TreeNode<T>* copyNode(const TreeNode<T>* node, TreeNode<T>* parent = nullptr) const {
    if (node == nullptr) {
      return nullptr;
    }

    // Create a new node with the same data
    TreeNode<T>* newNode = new TreeNode<T>(node->data);
    newNode->level = node->level;
    newNode->index = node->index;
    newNode->parent = parent;

    // Recursively copy the children
    for (const auto& child : node->children) {
      TreeNode<T>* copiedChild = copyNode(child, newNode);
      newNode->children.push_back(copiedChild);
    }

    return newNode;
  }

  void print(int x, int y, int ySpace, TreeNode<T>* node, std::vector<std::vector<int>> &array, int type);

public:
  int nNodes = 0;

  void print();

  Tree() : root(nullptr) {}

  // Copy constructor
  Tree(const Tree& original) {
    // Deep copy the entire tree structure
    if (original.root != nullptr) {
      root = copyNode(original.root);
    }
    else {
      root = nullptr;
    }
  }

  // Move constructor
  Tree(Tree&& other) noexcept {
    root = other.root; // Steal the root pointer
    nNodes = other.nNodes;
    other.root = nullptr; // Set the source object's root to nullptr
  }

  // Move assignment operator
  Tree& operator=(Tree&& other) noexcept {
    if (this != &other) { // Prevent self-assignment
      destroyTree(root); // Clean up current resources
      root = other.root; // Steal the source object's root
      other.root = nullptr; // Clear the source object's root
    }
    return *this;
  }

  // Copy assignment operator
  Tree& operator=(const Tree& original) {
    if (this != &original) {
      destroyTree(root); // Clean up current resources
      root = copyNode(original.root); // Deep copy the original tree
    }
    return *this;
  }

  ~Tree() {
    destroyTree(root);
  }

  TreeNode<T>* getRoot() const {
    return root;
  }

  void setRoot(TreeNode<T>* node) {
    root = node;
    node->index = 0;
    nNodes = 1;
  }

  void addChild(TreeNode<T>* parent, TreeNode<T>* child);

  void addChild(TreeNode<T>* parent, int nChildren);

  void getNodes(TreeNode<T>* node, std::vector<TreeNode<T>*> &allNodes) const {
    if(node != nullptr) {
      allNodes.push_back(node);
      for (TreeNode<T>* child : node->children) {
        getNodes(child, allNodes);
      }
    }
  }

  std::vector<TreeNode<T>*> getNodes() const {
    std::vector<TreeNode<T>*> allNodes;
    getNodes(root, allNodes);
    return allNodes;
  }

  std::vector<TreeNode<T>*> getLevel(int level) const {
    std::vector<TreeNode<T>*> allNodes = this->getNodes();
    std::vector<TreeNode<T>*> levelNodes;
    for(auto& node : allNodes) {
      if(node->level == level) {
        levelNodes.push_back(node);
      }
    }
    return levelNodes;
  }

  int depth() const {
    int output = 0;
    std::vector<TreeNode<T>*> level = getLevel(0);
    while(level.size() > 0) {
      output++;
      level = getLevel(output);
    }
    return output;
  }
};

#endif