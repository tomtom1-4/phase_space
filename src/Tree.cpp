#include "Tree.hpp"
#include "PhaseSpace.hpp"

template <typename T>
void Tree<T>::addChild(TreeNode<T>* parent, TreeNode<T>* child) {
  // If no root selected, the child will become the root
  if (root == nullptr) {
    root = child;
    parent = nullptr;
    return;
  }
  parent->children.push_back(child);
  if(child != nullptr) {
    child->parent=parent;
    child->level = parent->level + 1;
  }
}

template <typename T>
void Tree<T>::print(int x, int y, int ySpace, TreeNode<T>* node, std::vector<std::vector<int>> &array) {
  if(node == nullptr) return;
  array[x][y] = 1;
  for(int y_counter = y - ySpace/2; y_counter < y + ySpace/2; y_counter++) {
    array[x+1][y_counter] = 2;
  }
  for(int i = 0; i < node->children.size(); i++) {
    int new_ySpace = ySpace/node->children.size();
    print(x + 2, y + ySpace/2 - i*ySpace/(node->children.size() - 1), new_ySpace, node->children[i], array);
  }
}

template <typename T>
void Tree<T>::print() {
  int size = 50;
  std::vector<std::vector<int>> array(size, std::vector<int>(size, 0));
  if(root == nullptr) return;
  print(0, size/2, size/2, root, array);

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

/*template <typename T>
TreeNode<T> generate_TreeNode(TreeNode<T>& parent=nullptr, int nChildren) {
  if(nChildren == 2) {
    TreeNode<T>* child1 = new TreeNode<T>(T());
    TreeNode<T>* child2 = new TreeNode<T>(T());
    child1->parent = parent;
    child2->parent = parent;
    parent.children.push_back(child1);
    parent.children.push_back(child2);
    return parent;
  }

  for(int i = 2; i <= nChildren; i++) {
    return generate_TreeNode(parent, i);
  }
}

template <typename T>
std::vector<Tree<T>> generate_Tree(int nReference, int nUnresolved) {
  std::vector<Tree<T>> output;

  Tree<T> tree;
  TreeNode<T>* rootc = new TreeNode<T>(T());
  tree.setRoot(rootc);
  for(int i = 0; i < nReference; i++) {
    TreeNode<T>* referenceNode = new TreeNode<T>(T());
    tree.addChild(rootc, referenceNode);
  }
  for(int i = 0)
}*/
template <>
void Tree<Cluster>::print(int x, int y, int ySpace, TreeNode<Cluster>* node, std::vector<std::vector<int>> &array) {
  if(node == nullptr) return;
  array[x][y] = 1;
  for(int y_counter = y - ySpace/2; y_counter < y + (ySpace/2); y_counter++) {
    array[x+1][y_counter] = 2;
  }

  if(node->children.size() == 0) {
    array[x+2][y-ySpace/2] = 3;
    for(int i = 0; i < (node->data.unresolved); i++) {
      array[x+2][y + ySpace/2 - i*ySpace/(node->data.unresolved)] = 1;
    }
  }
  else {
    for(int i = 0; i < (node->children.size() + node->data.unresolved); i++) {
      int new_ySpace = ySpace/(node->children.size() + node->data.unresolved);
      if(i < node->children.size())
        print(x + 2, y + ySpace/2 - i*ySpace/((node->children.size() + node->data.unresolved) - 1), new_ySpace, node->children[i], array);
      else
        array[x+2][y + ySpace/2 - i*ySpace/((node->children.size() + node->data.unresolved) - 1)] = 3;
    }
  }
}

template <>
void Tree<Cluster>::print() {
  int size = 100;
  std::vector<std::vector<int>> array(50, std::vector<int>(size, 0));
  if(root == nullptr) return;
  print(0, size/2, size/2, root, array);

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

