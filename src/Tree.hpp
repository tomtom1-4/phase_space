#ifndef TREE_HPP
#define TREE_HPP

#include <iostream>
#include <vector>
#include <string>

template <typename T>
struct TreeNode {
  T data;
  int level = 0;
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
    newNode->parent = parent;

    // Recursively copy the children
    for (const auto& child : node->children) {
      TreeNode<T>* copiedChild = copyNode(child, newNode);
      newNode->children.push_back(copiedChild);
    }

    return newNode;
  }

  void print(int x, int y, int ySpace, TreeNode<T>* node, std::vector<std::vector<int>> &array);

public:
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

  ~Tree() {
    destroyTree(root);
  }

  TreeNode<T>* getRoot() const {
    return root;
  }

  void setRoot(TreeNode<T>* node) {
    root = node;
  }

  void addChild(TreeNode<T>* parent, TreeNode<T>* child);

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

  /*void preorderTraversal(TreeNode<T>* node) {
    if (node != nullptr) {
      if(node->parent != nullptr) {
        if((node->parent->children.size() > 1) and (node == node->parent->children[0])) std::cout << "{";
      }
      else {
        std::cout << "{";
      }
      std::cout << node->data;
      if(node->children.size() != 0) std::cout << ",";
      for (TreeNode<T>* child : node->children) {
        preorderTraversal(child);
      }
      if(node->parent != nullptr) {
        if((node->parent->children.size() > 1) and (node == node->parent->children.back())) std::cout << "}";
      }
      else {
        std::cout << "}";
      }
      if(node->parent != nullptr) {
        if(node != node->parent->children.back()) std::cout << ",";
      }
    }
  }


  void preorderTraversal() {
    preorderTraversal(root);
  }*/
};

//template <typename T>
//Tree<T> generate_Tree(int nReference, int nUnresolved);

#endif