#ifndef TREE_HPP
#define TREE_HPP

#include <iostream>
#include <vector>

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

public:
  Tree() : root(nullptr) {}

  ~Tree() {
    destroyTree(root);
  }

  TreeNode<T>* getRoot() const {
    return root;
  }

  void setRoot(TreeNode<T>* node) {
    root = node;
  }

  void addChild(TreeNode<T>* parent, TreeNode<T>* child) {
    // If no root selected, the child will become the root
    if (root == nullptr) {
      root = child;
      parent = nullptr;
      return;
    }
    parent->children.push_back(child);
    child->parent=parent;
    child->level = parent->level + 1;
  }

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

  void preorderTraversal(TreeNode<T>* node) {
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
  }
};

#endif