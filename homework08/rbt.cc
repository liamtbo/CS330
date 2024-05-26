#include "rbt.h"

// ---------------------------------------
// Node class
// Default constructor
RBTNode::RBTNode() : Node() {
    color = BLACK;
}

// Constructor
RBTNode::RBTNode(int in) : Node(in) {
    color = BLACK;
}
// Destructor
RBTNode::~RBTNode() {
}

void RBTNode::add_color(Node_color c)
{
  color = c;
}

void RBTNode::print_color(ostream& to)
{
    if(color == RED) {
        to << "Red";
    } else if (color == BLACK) {
        to << "Black";
    } else {
        cerr << "ERR: invalid color" << endl;
    }
}
void RBTNode::print_info(ostream& to)
{
    to << get_key() << " ";
    print_color(to);
}

Node_color RBTNode::get_color()
{
  return color;
}
// ---------------------------------------


// ---------------------------------------
// RBT class
// Constructor and destructor
RBT::RBT() : BST()
{
  sentinel = new RBTNode(-1);
  root = sentinel;
}
RBT::~RBT()
{
    // re-using BST's inorder_free
    inorder_free(root, sentinel);
    // This makes sure that root is set to nullptr, so that the parent class's
    // constructor does not try to free the tree again
    root = nullptr;
    delete sentinel; 
}

// Functions that are basically wrappers for the parent class functions
// minimum key in the BST
RBTNode* RBT::tree_min()
{
    // return (RBTNode*) get_min(root);
    // return (RBTNode*) BST::tree_min();
    return (RBTNode*) BST::get_min(root, sentinel);
}
// maximum key in the BST
RBTNode* RBT::tree_max()
{
    // return (RBTNode*) get_max(root);
    // return (RBTNode*) BST::tree_max();
    return (RBTNode*) BST::get_max(root, sentinel);
}

// Get successor of the given node
RBTNode* RBT::get_succ(RBTNode* in)
{
  return (RBTNode*) BST::get_succ((Node*) in, sentinel);
}

// Get predecessor of the given node
RBTNode* RBT::get_pred(RBTNode* in)
{
  return (RBTNode*) BST::get_pred((Node*) in, sentinel);
}

// Search the tree for a given key
RBTNode* RBT::tree_search(int search_key)
{
  return (RBTNode*) BST::tree_search(search_key, sentinel);
}

void RBT::walk(ostream& to)
{
  BST::inorder_walk(root, to, sentinel);
}



// New functions to RBT
// right rotate
void RBT::right_rotate(RBTNode* y)
{
    /* TODO */
    RBTNode *x = (RBTNode*) (y->get_left());

    y->add_left((RBTNode*) (x->get_right()));
    if ((RBTNode*) (x->get_right()) != this->get_sentinel()) {
        x->get_right()->add_parent(y);
    }
    x->add_parent((RBTNode*) (y->get_parent()));
    if ((RBTNode*) (y->get_parent()) == this->get_sentinel()) {
        this->root = x;
    }
    else if (y == (RBTNode*) (y->get_parent()->get_right())) {
        y->get_parent()->add_right(x);
    }
    else {
        y->get_parent()->add_left(x);
    }
    x->add_right(y);
    y->add_parent(x);
}
// Left rotate
void RBT::left_rotate(RBTNode* x)
{
    /* TODO */
    RBTNode *y = (RBTNode*) (x->get_right());
    x->add_right((RBTNode*) (y->get_left()));
    if ((RBTNode*) (y->get_left()) != this->get_sentinel()) {
        y->get_left()->add_parent(x);
    }
    y->add_parent((RBTNode*)(x->get_parent()));
    if ((RBTNode*) (x->get_parent()) == this->get_sentinel()) {
        this->root = y;
    }
    else if (x == (RBTNode*) (x->get_parent()->get_left())) {
        x->get_parent()->add_left(y);
    }
    else {
        x->get_parent()->add_right(y);
    }
    y->add_left(x);
    x->add_parent(y);
}
void RBT::rb_insert_fixup(RBTNode* in)
{
    /* TODO */
    RBTNode* x = in;
    while (x != root && ((RBTNode *)(x->get_parent()))->get_color() == RED) {
        if (x->get_parent() == x->get_parent()->get_parent()->get_left()) {
            RBTNode* y = (RBTNode*)(x->get_parent()->get_parent()->get_right());
            // if uncle and parent are red
            if (y->get_color() == RED) {
                ((RBTNode*)(x->get_parent()))->add_color(BLACK);
                y->add_color(BLACK);
                ((RBTNode*)(x->get_parent()->get_parent()))->add_color(RED);
                x = (RBTNode*)(x->get_parent()->get_parent());
            } // uncle is black and x is right child
            else {
                if (x == (RBTNode*) (x->get_parent()->get_right())) {
                    x = (RBTNode*)(x->get_parent());
                    left_rotate(x);
                }
                ((RBTNode*)(x->get_parent()))->add_color(BLACK);
                ((RBTNode*)(x->get_parent()->get_parent()))->add_color(RED);
                right_rotate(((RBTNode*)(x->get_parent()->get_parent())));
            }
        }
        else {
            RBTNode* y = (RBTNode*)(x->get_parent()->get_parent()->get_left());
            // if uncle and parent are red
            if (y->get_color() == RED) {
                ((RBTNode*)(x->get_parent()))->add_color(BLACK);
                y->add_color(BLACK);
                ((RBTNode*)(x->get_parent()->get_parent()))->add_color(RED);
                x = (RBTNode*)(x->get_parent()->get_parent());
            } // uncle is black and x is right child
            else {
                if (x == (RBTNode*) (x->get_parent()->get_left())) {
                    x = (RBTNode*)(x->get_parent());
                    right_rotate(x);
                }
                ((RBTNode*)(x->get_parent()))->add_color(BLACK);
                ((RBTNode*)(x->get_parent()->get_parent()))->add_color(RED);
                left_rotate(((RBTNode*)(x->get_parent()->get_parent())));
            }
        }
    }
    ((RBTNode*)(root))->add_color(BLACK);
}

void RBT::rb_insert_node(RBTNode* in)
{
    /* TODO */
    RBTNode* x = in;
    this->insert_node(x, sentinel);
    x->add_left(sentinel);
    x->add_right(sentinel);
    x->add_color(RED);
    if (x == root) {
        x->add_parent(sentinel);
    }
    rb_insert_fixup(x);
 
}

// ---------------------------------------

void RBT::rb_delete_fixup(RBTNode* in)
{
    /* TODO */
    RBTNode* x = in;
    while (x != root && x->get_color() == BLACK) {
        RBTNode* parent = (RBTNode*)x->get_parent();
        if (x == parent->get_left()) {
            RBTNode *w = (RBTNode*)(parent->get_right());
            // case 1
            if (w->get_color() == RED) {
                w->add_color(BLACK);
                parent->add_color(RED);
                left_rotate(parent);
                w = (RBTNode*)(parent->get_right());
            }
            // case 2
            if (((RBTNode*)(w->get_left()))->get_color() == BLACK \
                && ((RBTNode*)(w->get_right()))->get_color() == BLACK) { // here
                w->add_color(RED);
                x = parent;
            }
            // case 3-4
            else {
                if (((RBTNode*)(w->get_right()))->get_color() == BLACK) {
                    ((RBTNode*)(w->get_left()))->add_color(BLACK);
                    w->add_color(RED);
                    right_rotate(w);
                    w = (RBTNode*)(parent->get_right());
                }
                w->add_color(parent->get_color());
                parent->add_color(BLACK);
                ((RBTNode*)(w->get_right()))->add_color(BLACK);
                left_rotate(parent);
                x = (RBTNode*)root;
            }
        }
        else {
            RBTNode *w = (RBTNode*)(parent->get_left());
            // case 1
            if (w->get_color() == RED) {
                w->add_color(BLACK);
                parent->add_color(RED);
                right_rotate(parent);
                w = (RBTNode*)(parent->get_left());
            }
            // case 2
            if (((RBTNode*)(w->get_right()))->get_color() == BLACK \
                && ((RBTNode*)(w->get_left()))->get_color() == BLACK) { // here
                w->add_color(RED);
                x = parent;
            }
            // case 3-4
            else {
                if (((RBTNode*)(w->get_left()))->get_color() == BLACK) {
                    ((RBTNode*)(w->get_right()))->add_color(BLACK);
                    w->add_color(RED);
                    left_rotate(w);
                    w = (RBTNode*)(parent->get_left());
                }
                w->add_color(parent->get_color());
                parent->add_color(BLACK);
                ((RBTNode*)(w->get_left()))->add_color(BLACK);
                right_rotate(parent);
                x = (RBTNode*)root;
            }
        }
    }
    x->add_color(BLACK);
}

void RBT::rb_delete_node(RBTNode* out)
{
    /* TODO */
    // temporary pointers
    RBTNode* y;
    RBTNode* x;

    // first determine if out has at most one child
    if(((RBTNode*)(out->get_left()) == sentinel) || ((RBTNode*)(out->get_right()) == sentinel)) {
        // 0 or 1 child - out will be spliced
        y = out;
    } else {
        // 2 children - out's successor will be sliced out
        y = this->get_succ(out);
        // cout << y->get_key() << endl;
    }
    Node_color orig_color = y->get_color();
    // y = out -> 1 child, succ can't have left child
    RBTNode *parent = (RBTNode*)(y->get_parent());
    if((RBTNode*)(y->get_left()) != sentinel) {
        x = (RBTNode*)(y->get_left());                            
    } // y = out or y = succ(y) -> both only have the right child
    else { 
        x = (RBTNode*)(y->get_right());
    }

    // Splice out node y with its child
    // check for null, otherwise the code will segfault
    // if(x != sentinel) {
    //     x->add_parent((RBTNode*)(y->get_parent()));
    // }

    // Extra code to account for RBT
    if(x == sentinel && x != nullptr) { // (perror == -> !=)
        x->add_parent(parent);
    }

    if(parent == sentinel) {
        root = x;
    } else {
        if(y == parent->get_left()) {
            parent->add_left(x);
        } else {
            parent->add_right(x);
        }
    }
    // Finally, replace out with y
    if(y != out) {
        // first copy out's parent and children to y
        y->add_parent(out->get_parent());
        if(y->get_parent() == sentinel) {
            root = y;
        } else {
            // update the parent's child to y
            if(out == (out->get_parent()->get_left())) {
                y->get_parent()->add_left(y);
            } else {
                y->get_parent()->add_right(y);
            }
        }
        y->add_left((RBTNode*)(out->get_left()));
        y->add_right((RBTNode*)(out->get_right()));
        // y->add_color(((RBTNode*)(out->get_left()))->get_color());
        // y->add_color(orig_color);

        // update the children's parents to y
        y->get_left()->add_parent(y);
        y->get_right()->add_parent(y);
        y->add_color(out->get_color());
    }

    delete out;
    if (orig_color == BLACK) {
        rb_delete_fixup(x);
    }
}
// ---------------------------------------

