#include "bst.h"

// ---------------------------------------
// Node class
// Default constructor
Node::Node() {
    // TODO: Implement this
    key = 0;
    parent = nullptr;
    right = nullptr;
    left = nullptr;
}
// Constructor
Node::Node(int in) {
    // TODO: Implement this
    key = in;
    parent = nullptr;
    right = nullptr;
    left = nullptr;
}
// Destructor
Node::~Node() {
    // TODO: Implement this
}

// Add parent 
void Node::add_parent(Node* in) {
    // TODO: Implement this
    this->parent = in;
}
// Add to left of current node
void Node::add_left(Node* in) {
    // TODO: Implement this
    this->left = in;
}
// Add to right of current node
void Node::add_right(Node* in) {
    // TODO: Implement this
    this->right = in;
}

// Get key
int Node::get_key()
{
    // TODO: Implement this
    return this->key;
}
// Get parent node
Node* Node::get_parent()
{
    // TODO: Implement this
    return this->parent;
}
// Get left node
Node* Node::get_left()
{
    // TODO: Implement this
    return this->left;
}
// Get right node
Node* Node::get_right()
{
    // TODO: Implement this
    return this->right;
}
// Print the key to ostream to
// Do not change this
void Node::print_info(ostream& to)
{
    to << key << endl;
}
// ---------------------------------------


// ---------------------------------------
// BST class
// Walk the subtree from the given node
void BST::inorder_walk(Node* in, ostream& to)
{
    // TODO: Implement this
    if (in) {
        this->inorder_walk(in->get_left(), to);
        in->print_info(to);
        this->inorder_walk(in->get_right(), to);
    }
}
// Constructor
BST::BST()
{
    // TODO: Implement this
    root = nullptr;
}
// Destructor
BST::~BST()
{
    // TODO: Implement this
    while (root != nullptr) {
        this->delete_node(this->tree_min());
    }
}
// Insert a node to the subtree
void BST::insert_node(Node* in)
{
    // TODO: Implement this
   
    Node* y = nullptr;
    Node* x = this->root;
    while (x != nullptr) {
        // saves parent
        y = x;
        // tests for when new leaf is hit
        if (in->get_key() < x->get_key()) {
            x = x->get_left();
        } else {
            x = x->get_right();
        }
    }
    in->add_parent(y);
    if (y == nullptr) {
        this->root = in;
    } else if (in->get_key() < y->get_key()) {
        y->add_left(in);
    } else {
        y->add_right(in);
    }
}
// Delete a node to the subtree
void BST::delete_node(Node* out)
{
    // TODO: Implement this
    Node *x;
    Node *y;
    Node *z = out;
    if (z->get_left() == nullptr || z->get_right() == nullptr) {
        y = z;
    } else {
        y = this->get_succ(z);
    }
    if (y->get_left() != nullptr) {
        x = y->get_left();
    } else {
        x = y->get_right();
    }
    if (x != nullptr) {
        x->add_parent(y->get_parent());
    }
    if (y->get_parent() == nullptr) {
        this->root = x;
    } else if (y == y->get_parent()->get_left()) {
        y->get_parent()->add_left(x);
    } else {
        y->get_parent()->add_right(x);
    }
    // need to copy key over but there's no built in function
    if (y != z) {
        Node cpy(y->get_key());
        cpy.add_left(z->get_left());
        cpy.add_right(z->get_right());
        cpy.add_parent(z->get_parent());
        *z = cpy;
    }
    delete y;
}
// minimum key in the BST
Node* BST::tree_min()
{
    // TODO: Implement this
    Node* x = root;
    while(x->get_left() != nullptr) {
        x = x->get_left();
    }
    return x;
}
// maximum key in the BST
Node* BST::tree_max()
{
    // TODO: Implement this
    Node* x = root;
    while (x->get_right() != nullptr) {
        x = x->get_right();
    }
    return x;
}
// Get the minimum node from the subtree of given node
Node* BST::get_min(Node* in)
{
    // TODO: Implement this
    Node* x = in;
    while(x->get_left() != nullptr) {
        x = x->get_left();
    }
    return x;
}
// Get the maximum node from the subtree of given node
Node* BST::get_max(Node* in)
{
// TODO: Implement this
    Node* x = in;
    while (x->get_right() != nullptr) {
        x = x->get_right();
    }
    return x;
}
// Get successor of the given node
Node* BST::get_succ(Node* in)
{
    // TODO: Implement this
    // Node with the next largest key
    Node *x = in;
    if (x->get_right() != nullptr) {
        return get_min(x->get_right());
    }
    Node* y = x->get_parent();
    while (y != nullptr && x == y->get_right()) {
        x = y;
        y = y->get_parent();
    }
    return y;
}
// Get predecessor of the given node
Node* BST::get_pred(Node* in)
{
    // TODO: Implement this
    // Node with the next smallest key
    Node *x = in;
    if (x->get_left() != nullptr) {
        return get_max(x->get_left());
    }
    Node *y = x->get_parent();
    while (y != nullptr && x == y->get_left()) {
        x = y;
        y = y->get_parent();
    }
    return y;
}
// Walk the BST from min to max
void BST::walk(ostream& to)
{
    // TODO: Implement this
    inorder_walk(root, to);
}
// Search the tree for a given key
Node* BST::tree_search(int search_key)
{
    // TODO: Implement this
    Node *x = root;
    while (x != nullptr && search_key != x->get_key()) {
        if (search_key < x->get_key()) {
            x = x->get_left();
        } else {
            x = x->get_right();
        }
    }
    return x;
}
// ---------------------------------------
