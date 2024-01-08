#ifndef SRC_SORTERNODE_H_
#define SRC_SORTERNODE_H_
#include "util.h"
template<class T>
class sorter_node {
public:
	static sorter_node* insert(sorter_node* root, const T& val, int from) {
		if (root) {
			bool tmp;
			return root->_M_insert(val, from, tmp);
		} else {
			sorter_node* rtn = (sorter_node*)malloc(sizeof(sorter_node));
			rtn->_M_val = val;
			rtn->_M_from = from;
			rtn->_M_bal = 0;
			rtn->_M_l_node = rtn->_M_r_node = NULL;
			return rtn;
		}
	}
/*
	sorter_node* insert(T val, int from) {
		bool tmp;
		return _M_insert(val, from, tmp);
	}
*/
	sorter_node* get_min(T *val, int &from) {
		bool tmp;
		return _M_get_min(val, from, tmp);
	}
private:
	T _M_val;
	int _M_from;
	char _M_bal;
	sorter_node *_M_l_node, *_M_r_node;
	sorter_node* _M_insert(T val, int from, bool &changed) {
		if (_M_val < val) {
			if (!_M_r_node) {
				_M_r_node = sorter_node::insert(NULL, val, from);
				changed = _M_bal++ == 0;
				return this;
			}
			_M_r_node = _M_r_node->_M_insert(val, from, changed);
			if (!changed) return this;
			if (_M_bal++ <= 0) {
				changed = _M_bal == 1;
				return this;
			}
			if (_M_r_node->_M_bal > 0) {
				sorter_node *root = _M_r_node;
				_M_r_node = root->_M_l_node;
				root->_M_l_node = this;
				root->_M_bal = _M_bal = 0;
				changed = false;
				return root;
			} else {
				sorter_node *root = _M_r_node->_M_l_node;
				sorter_node *root_r_node = root->_M_r_node;
				sorter_node *root_l_node = root->_M_l_node;
				root->_M_r_node = _M_r_node;
				root->_M_l_node = this;
				root->_M_r_node->_M_l_node = root_r_node;
				_M_r_node = root_l_node;
				if (root->_M_bal == 0) 
					root->_M_l_node->_M_bal = root->_M_r_node->_M_bal = 0;
				else if (root->_M_bal == 1) {
					root->_M_l_node->_M_bal = -1;
					root->_M_r_node->_M_bal = 0;
				} else {
					root->_M_l_node->_M_bal = 0;
					root->_M_r_node->_M_bal = 1;
				}
				root->_M_bal = 0;
				changed = false;
				return root;
			}
		} else {
			if (!_M_l_node) {
				_M_l_node = sorter_node::insert(NULL, val, from);
				changed = _M_bal-- == 0;
				return this;
			}
			_M_l_node = _M_l_node->_M_insert(val, from, changed);
			if (!changed) return this;
			if (_M_bal-- >= 0) {
				changed = _M_bal == -1;
				return this;
			}
			if (_M_l_node->_M_bal < 0) {
				sorter_node *root = _M_l_node;
				_M_l_node = root->_M_r_node;
				root->_M_r_node = this;
				root->_M_bal = _M_bal = 0;
				changed = false;
				return root;
			} else {
				sorter_node *root = _M_l_node->_M_r_node;
				sorter_node *root_l_node = root->_M_l_node;
				sorter_node *root_r_node = root->_M_r_node;
				root->_M_l_node = _M_l_node;
				root->_M_r_node = this;
				root->_M_l_node->_M_r_node = root_l_node;
				_M_l_node = root_r_node;
				if (root->_M_bal == 0) 
					root->_M_l_node->_M_bal = root->_M_r_node->_M_bal = 0;
				else if (root->_M_bal == 1) {
					root->_M_l_node->_M_bal = -1;
					root->_M_r_node->_M_bal = 0;
				} else {
					root->_M_l_node->_M_bal = 0;
					root->_M_r_node->_M_bal = 1;
				}
				root->_M_bal = 0;
				changed = false;
				return root;
			}
		}
	}
	sorter_node* _M_get_min(T *val, int &from, bool &changed) {
		if (!_M_l_node) {
			*val = _M_val;
			from = _M_from;
			sorter_node *rtn = _M_r_node;
			delete this;
			changed = true;
			return rtn;
		} else {
			_M_l_node = _M_l_node->_M_get_min(val, from, changed);
			if (!changed) return this;
			if (_M_bal++ <= 0) {
				changed = _M_bal != 1;
				return this;
			}
			if (_M_r_node->_M_bal >= 0) {
				sorter_node *root = _M_r_node;
				_M_r_node = root->_M_l_node;
				root->_M_l_node= this;
				changed = root->_M_bal > 0; 
				if (root->_M_bal > 0) 
					root->_M_bal = _M_bal = 0;
				else {
					root->_M_bal = -1;
					_M_bal = 1;
				}
				return root;
			} else {
				sorter_node *root = _M_r_node->_M_l_node;
				sorter_node *root_r_node = root->_M_r_node;
				sorter_node *root_l_node = root->_M_l_node;
				root->_M_r_node = _M_r_node;
				root->_M_l_node = this;
				root->_M_r_node->_M_l_node = root_r_node;
				_M_r_node = root_l_node;
				if (root->_M_bal == 0) 
					root->_M_l_node->_M_bal = root->_M_r_node->_M_bal = 0;
				else if (root->_M_bal == 1) {
					root->_M_l_node->_M_bal = -1;
					root->_M_r_node->_M_bal = 0;
				} else {
					root->_M_l_node->_M_bal = 0;
					root->_M_r_node->_M_bal = 1;
				}
				root->_M_bal = 0;
				changed = true;
				return root;
			}
		}
	}
};
#endif
