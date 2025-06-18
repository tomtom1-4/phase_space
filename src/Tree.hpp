#ifndef TREE_HPP
#define TREE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace PSF
{

/**
  * @struct TreeNode
  * @brief A node in a generic tree structure, storing data of type T.
  *
  * Each node holds a value of type T and pointers to its parent and children.
  * Children are stored as std::unique_ptr<TreeNode<T>> to automate memory management.
  *
  * @tparam T The data type stored in each node.
  */
template <typename T>
struct TreeNode {
  /**
    * @brief The data stored in this node.
    */
  T data;

  /**
    * @brief The depth level of this node in the tree:
    *        0 for the root, 1 for its children, etc.
    */
  int level = 0;

  /**
    * @brief A unique index for this node among its siblings, starting at 0.
    */
  int index = 0;

  /**
    * @brief Pointers to this node's children as unique_ptr.
    *
    * When this node is destroyed, all children are also destroyed automatically.
    */
  std::vector<std::unique_ptr<TreeNode<T>>> children;

  /**
    * @brief Pointer to this node's parent. This is a raw pointer because the parent
    *        does not own the child. It may be nullptr if this node is the root.
    */
  TreeNode<T>* parent = nullptr;

  /**
    * @brief Constructs a TreeNode with the specified initial data.
    * @param value The data to store in the new node.
    */
  explicit TreeNode(const T& value)
    : data(value)
  {}
};

/**
  * @class Tree
  * @brief A generic tree data structure with smart-pointer-based ownership.
  *
  * This class provides a root node and various utility functions for adding children,
  * traversing levels, getting all nodes, copying, and printing. It supports deep copy
  * and move semantics. Memory is automatically cleaned up when the Tree is destroyed.
  *
  * @tparam T The data type of each node in the tree.
  */
template <typename T>
class Tree {
private:
  /**
    * @brief The root node of the tree, or nullptr if empty.
    *
    * Owned by the Tree via a unique_ptr; when this pointer is reset or destroyed,
    * the entire tree is freed.
    */
  std::unique_ptr<TreeNode<T>> m_root;

  /**
    * @brief The total number of nodes in the tree.
    *
    * Set to 0 when the tree is empty. Updated whenever
    * a node is added or removed.
    */
  int m_nodeCount = 0;

  /**
    * @brief Recursively copies a node and all its descendants.
    *
    * @param node The node to copy from (source).
    * @param parent The new parent node in the copy, or nullptr if this is the new root.
    * @return A unique_ptr to the newly created TreeNode<T>.
    */
  std::unique_ptr<TreeNode<T>> copyNode(const TreeNode<T>* node,
                                        TreeNode<T>* parent = nullptr) const
  {
    if (!node) {
      return nullptr;
    }

    // Create a new node with the same data
    std::unique_ptr<TreeNode<T>> newNode = std::make_unique<TreeNode<T>>(node->data);
    newNode->level = node->level;
    newNode->index = node->index;
    newNode->parent = parent;

    // Recursively copy each child
    newNode->children.reserve(node->children.size());
    for (const auto& child : node->children) {
      std::unique_ptr<TreeNode<T>> copiedChild = copyNode(child.get(), newNode.get());
      newNode->children.push_back(std::move(copiedChild));
    }
    return newNode;
  }

  /**
    * @brief Recursive helper for printing the tree (skeleton function).
    *
    * @param x A horizontal position parameter.
    * @param y A vertical position parameter.
    * @param ySpace The vertical spacing used in printing.
    * @param node The current node in the traversal.
    * @param array A 2D buffer or structure for representing the tree visually.
    * @param type A formatting parameter or style indicator.
    *
    * @note You would implement the details in this function as needed, or
    *       move it to a .cpp/.ipp file. For now, it is just declared.
    */
  void print(int x, int y, int ySpace, TreeNode<T>* node,
             std::vector<std::vector<int>>& array, int type);

  /**
    * @brief Recursively traverses the subtree to gather all nodes.
    * @param node The current node being visited.
    * @param allNodes A reference to a vector collecting the visited nodes.
    */
  void getNodesRec(TreeNode<T>* node, std::vector<TreeNode<T>*>& allNodes) const
  {
    if (node) {
      allNodes.push_back(node);
      for (const auto& child : node->children) {
        getNodesRec(child.get(), allNodes);
      }
    }
  }

public:
  /**
    * @brief Default constructor, creating an empty tree.
    */
  Tree() = default;

  /**
    * @brief Copy constructor that performs a deep copy of an existing tree.
    * @param original The tree to copy from.
    */
  Tree(const Tree& original)
  {
    if (original.m_root) {
      m_root = copyNode(original.m_root.get(), nullptr);
    }
    m_nodeCount = original.m_nodeCount;
  }

  /**
    * @brief Move constructor that takes ownership of the source tree's data.
    * @param other The tree to move from.
    * @note After this operation, @p other becomes empty (root set to nullptr).
    */
  Tree(Tree&& other) noexcept
    : m_root(std::move(other.m_root)),
      m_nodeCount(other.m_nodeCount)
  {
    other.m_nodeCount = 0;
  }

  /**
    * @brief Copy assignment operator that deep copies the given tree.
    * @param original The tree to copy from.
    * @return A reference to this tree.
    */
  Tree& operator=(const Tree& original)
  {
    if (this != &original) {
      m_root.reset();
      if (original.m_root) {
        m_root = copyNode(original.m_root.get(), nullptr);
      }
      m_nodeCount = original.m_nodeCount;
    }
    return *this;
  }

  /**
    * @brief Move assignment operator that takes ownership of the source tree's data.
    * @param other The tree to move from.
    * @return A reference to this tree.
    */
  Tree& operator=(Tree&& other) noexcept
  {
    if (this != &other) {
      m_root.reset();
      m_root = std::move(other.m_root);
      m_nodeCount = other.m_nodeCount;
      other.m_nodeCount = 0;
    }
    return *this;
  }

  /**
    * @brief Destructor, automatically cleans up the entire tree.
    *
    * Since children are stored as std::unique_ptr, no manual destruction is required.
    */
  ~Tree() = default;

  /**
    * @brief Retrieves a pointer to the root node of this tree.
    * @return A raw pointer to the root node, or nullptr if the tree is empty.
    *
    * @warning The ownership of the root remains with this Tree, so don't delete it.
    */
  TreeNode<T>* getRoot() const
  {
    return m_root.get();
  }

  /**
    * @brief Dynamically creates and sets a new root node with the specified data.
    * @param value The data to store in the new root node.
    *
    * Resets any existing tree structure, since the old root is replaced.
    */
  void setRoot(const T& value)
  {
    m_root = std::make_unique<TreeNode<T>>(value);
    m_root->index = 0;
    m_root->level = 0;
    m_root->parent = nullptr;
    m_nodeCount = 1;
  }

  /**
    * @brief Sets the root node from an existing std::unique_ptr.
    * @param newRoot A unique_ptr to a node that becomes the new root.
    * @warning This replaces any existing root and tree structure.
    */
  void setRoot(std::unique_ptr<TreeNode<T>> newRoot)
  {
    m_root = std::move(newRoot);
    if (m_root) {
      m_root->index = 0;
      m_root->level = 0;
      m_root->parent = nullptr;
      m_nodeCount = 1; // Adjust as needed if you build a larger subtree before calling this.
    } else {
      m_nodeCount = 0;
    }
  }

  /**
    * @brief Creates and adds a new child node under a given parent, initializing it with data.
    * @param parent A raw pointer to the parent node (which must belong to this tree).
    * @param child The unique_ptr to the child to add to the parents children.
    *
    * Increments the node count. Sets the child's parent, level, and index.
    * Throws std::invalid_argument if parent is null.
    */
  void addChild(TreeNode<T>* parent, std::unique_ptr<TreeNode<T>> child);

  /**
    * @brief Creates and adds multiple new children under a given parent, each holding default-constructed data.
    * @param parent A pointer to the parent node. Must not be null.
    * @param nChildren The number of child nodes to add.
    *
    * Each new child has the same level and indexing is based on current size of the parent's children vector.
    */
  void addChild(TreeNode<T>* parent, int nChildren);

  /**
    * @brief Gathers pointers to all nodes in the tree via a pre-order traversal.
    * @return A std::vector of raw pointers to each node. May be empty if the tree is empty.
    *
    * @warning The ownership remains with the Tree. Do not store these pointers long-term without caution.
    */
  std::vector<TreeNode<T>*> getNodes() const
  {
    std::vector<TreeNode<T>*> allNodes;
    getNodesRec(m_root.get(), allNodes);
    return allNodes;
  }

  /**
    * @brief Retrieves all nodes at a particular depth level.
    * @param level The tree level to select (0 = root, 1 = children of root, etc.).
    * @return A std::vector of pointers to nodes at the given level.
    */
  std::vector<TreeNode<T>*> getLevel(int level) const
  {
    std::vector<TreeNode<T>*> allNodes = getNodes();
    std::vector<TreeNode<T>*> levelNodes;
    for (auto* node : allNodes) {
      if (node->level == level) {
        levelNodes.push_back(node);
      }
    }
    return levelNodes;
  }

  /**
    * @brief Computes the depth of the tree by counting how far levels extend.
    * @return The number of levels in the tree. Returns 0 if the tree is empty.
    */
  int depth() const
  {
    int d = 0;
    // Level 0 is root; keep checking levels until empty
    auto currentLevel = getLevel(d);
    while (!currentLevel.empty()) {
      d++;
      currentLevel = getLevel(d);
    }
    return d;
  }

  /**
    * @brief Returns the total number of nodes in this tree.
    * @return An integer count of nodes, or 0 if the tree is empty.
    */
  int size() const
  {
    return m_nodeCount;
  }

  /**
    * @brief Prints the tree.
    */
  void print();
};

}
#endif
