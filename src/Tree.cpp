#include "Tree.hpp"
#include "PhaseSpace.hpp"

namespace PSF
{

template <typename T>
void Tree<T>::addChild(TreeNode<T>* parent, TreeNode<T>* child) {
  // If no root selected, the child will become the root
  if (root == nullptr) {
    root = child;
    parent = nullptr;
    child->index = 0;
    nNodes = 1;
    return;
  }
  parent->children.push_back(child);
  if(child != nullptr) {
    child->parent=parent;
    child->level = parent->level + 1;
    child->index = nNodes;
    nNodes++;
  }
}

template <typename T>
void Tree<T>::addChild(TreeNode<T>* parent, int nChildren) {
  for(int i = 0; i < nChildren; i++) {
    TreeNode<T>* child = new TreeNode<T>(T());
    this->addChild(parent, child);
  }
}

template <typename T>
void Tree<T>::print(int x, int y, int ySpace, TreeNode<T>* node, std::vector<std::vector<int>> &array, int type) {
  if(node == nullptr) return;
  array[x][y] = 1;
  for(int y_counter = y - ySpace/2; y_counter < y + ySpace/2; y_counter++) {
    array[x+1][y_counter] = 2;
  }
  for(int i = 0; i < node->children.size(); i++) {
    int new_ySpace = ySpace/node->children.size();
    print(x + 2, y + ySpace/2 - i*ySpace/(node->children.size() - 1), new_ySpace, node->children[i], array, 0);
  }
}

template <typename T>
void Tree<T>::print() {
  int size = 50;
  std::vector<std::vector<int>> array(size, std::vector<int>(size, 0));
  if(root == nullptr) return;
  print(0, size/2, size/2, root, array, 1);

  for(int y = 0; y < size; y++) {
    for(int x = 0; x < size; x++) {
      switch (array[x][y])
      {
      case 0:
        std::cout << "   ";
        break;
      case 1:
        std::cout << "---";
        break;
      case 2:
        std::cout << " | ";
        break;
      default:
        break;
      }
    }
    std::cout << std::endl;
  }
}

template <>
void Tree<Cluster>::print(int x, int y, int ySpace, TreeNode<Cluster>* node, std::vector<std::vector<int>> &array, int type) {
  if(node == nullptr) return;
  array[x][y] = type;


  if(node->children.size() == 0) {
    return;
  }
  else {
    if(node->children.size() == 1) {
      print(x + 2, y, ySpace, node->children[0], array, 3);
    }
    else {
      for(int y_counter = y - ySpace/2; y_counter < y + (ySpace/2); y_counter++) {
        array[x+1][y_counter] = 2;
      }

      for(int i = 0; i < node->children.size(); i++) {
        int new_ySpace = ySpace/node->children.size();
        //if(i < node->data.unresolved)
        if(!node->children[i]->data.isReference)
          print(x + 2, y + ySpace/2 - i*ySpace/(node->children.size() - 1), new_ySpace, node->children[i], array, 1);
        else
          print(x + 2, y + ySpace/2 - i*ySpace/(node->children.size() - 1), new_ySpace, node->children[i], array, 3);
      }
    }
  }
}

template <>
void Tree<Cluster>::print() {
  int size = 100;
  std::vector<std::vector<int>> array(50, std::vector<int>(size, 0));
  if(root == nullptr) return;
  print(0, size/2, size/2, root, array, 4);

  for(int x = 0; x < array.size(); x++) {
    bool empty = true;
    for(int j = 0; j < array[0].size(); j++) {
      if(array[x][j] != 0) empty = false;
    }
    if(!empty) {
      for(int y = 0; y < array[0].size(); y++) {
        switch (array[x][y])
        {
        case 0:
          std::cout << " ";
          break;
        case 1:
          std::cout << "|"; // Unresolved
          break;
        case 2:
          std::cout << "-";
          break;
        case 3:
          std::cout << "#"; // Reference
          break;
        case 4:
          std::cout << "O"; // Root
        default:
          break;
        }
        //std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

template class Tree<Cluster>;

}
