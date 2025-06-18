#include "Tree.hpp"
#include "PhaseSpace.hpp"

namespace PSF
{

/**
  * @brief Creates and adds a new child node under a specified parent node.
  *
  * If the tree is empty (no root), the provided child becomes the new root.
  * Otherwise, the child is added to the parent's child list.
  *
  * @param parent A pointer to the parent node (non-owning).
  * @param child A std::unique_ptr to the newly created child. Ownership is transferred to the tree.
  * @throws std::invalid_argument if there is already a root but the @p parent is null.
  */
template <typename T>
void Tree<T>::addChild(TreeNode<T>* parent, std::unique_ptr<TreeNode<T>> child)
{
  // If there's no root, the child becomes the new root
  if (!m_root) {
    m_root = std::move(child);
    m_root->parent = nullptr;
    m_root->index = 0;
    m_root->level = 0;
    m_nodeCount = 1;
    return;
  }

  // Otherwise, attach the child under the existing parent
  if (!parent) {
    throw std::invalid_argument("Parent pointer is null, but the tree already has a root.");
  }

  // Move the unique_ptr into the parent's children vector
  parent->children.push_back(std::move(child));

  // Update the child's metadata
  TreeNode<T>* rawChild = parent->children.back().get();
  if (rawChild) {
    rawChild->parent = parent;
    rawChild->level = parent->level + 1;
    rawChild->index = m_nodeCount;
    m_nodeCount++;
  }
}

/**
  * @brief Creates and adds a specified number of default-constructed children under a given parent.
  *
  * Each child is constructed by calling T() (the default constructor for T).
  *
  * @param parent A raw pointer to the parent node.
  * @param nChildren The number of child nodes to create and add.
  */
template <typename T>
void Tree<T>::addChild(TreeNode<T>* parent, int nChildren)
{
  for (int i = 0; i < nChildren; i++) {
    // Create the child using a unique_ptr
    auto childPtr = std::make_unique<TreeNode<T>>(T());
    // Delegate to the unique_ptr overload
    addChild(parent, std::move(childPtr));
  }
}

/**
  * @brief Recursive helper function for 2D ASCII printing of the tree.
  *
  * @param x Current horizontal index in the 2D buffer.
  * @param y Current vertical index in the 2D buffer.
  * @param ySpace Vertical spacing used to position children.
  * @param node The current node being processed (non-owning pointer).
  * @param array A 2D buffer to mark lines and nodes.
  * @param type An integer that controls the drawn symbol (e.g., line or node).
  */
template <typename T>
void Tree<T>::print(int x, int y, int ySpace,
                    TreeNode<T>* node,
                    std::vector<std::vector<int>>& array,
                    int type)
{
  if (!node) {
    return;
  }
  // Mark this node's position in the array
  array[x][y] = 1;

  // Draw a vertical line segment
  for (int y_counter = y - ySpace / 2; y_counter < y + ySpace / 2; y_counter++) {
    array[x + 1][y_counter] = 2;
  }

  // Recursively print children
  for (int i = 0; i < static_cast<int>(node->children.size()); i++) {
    int new_ySpace = (ySpace == 0 || node->children.size() == 1)
                     ? ySpace
                     : ySpace / static_cast<int>(node->children.size());

    TreeNode<T>* childPtr = node->children[i].get();
    if (!childPtr) {
      continue;
    }

    int offset = (node->children.size() > 1)
                 ? i * (ySpace / (node->children.size() - 1))
                 : 0;

    print(x + 2,
          y + ySpace / 2 - offset,
          new_ySpace,
          childPtr,
          array,
          0);
  }
}

/**
  * @brief Prints the tree in ASCII form using a 2D grid (size=50).
  */
template <typename T>
void Tree<T>::print()
{
  int size = 50;
  std::vector<std::vector<int>> array(size, std::vector<int>(size, 0));

  // If empty, just return
  if (!m_root) {
    return;
  }

  // Begin printing from the middle
  print(0, size / 2, size / 2, m_root.get(), array, 1);

  // Render
  for (int y = 0; y < size; y++) {
    for (int x = 0; x < size; x++) {
      switch (array[x][y]) {
        case 0:
          std::cout << "   "; // blank
          break;
        case 1:
          std::cout << "---"; // node
          break;
        case 2:
          std::cout << " | "; // vertical line
          break;
        default:
          break;
      }
    }
    std::cout << std::endl;
  }
}

/**
  * @brief Specialization of print(int, int, int, TreeNode<Cluster>*, ...) for a Cluster tree.
  *
  * Distinguishes reference vs. unresolved clusters.
  * 'type' controls the ASCII symbol to display.
  */
template <>
void Tree<Cluster>::print(int x, int y, int ySpace,
                          TreeNode<Cluster>* node,
                          std::vector<std::vector<int>>& array,
                          int type)
{
  if (!node) {
    return;
  }
  array[x][y] = type;

  // Leaf node
  if (node->children.empty()) {
    return;
  }
  // Single child
  else if (node->children.size() == 1) {
    TreeNode<Cluster>* childPtr = node->children[0].get();
    if (childPtr) {
      print(x + 2, y, ySpace, childPtr, array, 3);
    }
  }
  // Multiple children
  else {
    for (int y_counter = y - ySpace / 2; y_counter < y + (ySpace / 2); y_counter++) {
      array[x + 1][y_counter] = 2;
    }

    for (int i = 0; i < static_cast<int>(node->children.size()); i++) {
      int new_ySpace = (ySpace == 0 || node->children.size() == 1)
                       ? ySpace
                       : (ySpace / static_cast<int>(node->children.size()));

      TreeNode<Cluster>* childPtr = node->children[i].get();
      if (!childPtr) {
        continue;
      }

      int offset = (node->children.size() > 1)
                   ? i * (ySpace / (node->children.size() - 1))
                   : 0;

      // '#'(3) for reference, '|'(1) for unresolved
      int symbol = (childPtr->data.isReference) ? 3 : 1;

      print(x + 2,
            y + ySpace / 2 - offset,
            new_ySpace,
            childPtr,
            array,
            symbol);
    }
  }
}

/**
  * @brief Specialized print() for Tree<Cluster>, using a larger 2D grid
  * and distinct ASCII symbols for reference vs. unresolved.
  */
template <>
void Tree<Cluster>::print()
{
  int width = 50;
  int height = 100;
  std::vector<std::vector<int>> array(width, std::vector<int>(height, 0));

  if (!m_root) {
    return;
  }

  // '4' used as an ASCII code for root node
  print(0, height / 2, height / 2, m_root.get(), array, 4);

  // Render each column that isn't empty
  for (int x = 0; x < width; x++) {
    bool empty = true;
    for (int y = 0; y < height; y++) {
      if (array[x][y] != 0) {
        empty = false;
        break;
      }
    }
    if (!empty) {
      for (int y = 0; y < height; y++) {
        switch (array[x][y]) {
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
            break;
          default:
            break;
        }
      }
      std::cout << std::endl;
    }
  }
}

// Explicitly instantiate the template class for Cluster so the linker generates code for it
template class Tree<Cluster>;

}