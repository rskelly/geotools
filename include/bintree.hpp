/*
 * bintree.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: rob
 */

#ifndef INCLUDE_BINTREE_HPP_
#define INCLUDE_BINTREE_HPP_

template <class K, class V>
class BinTree {
private:
	BinTree<K, V>* left;
	BinTree<K, V>* right;

	inline K abs(K a) {
		return a > 0 ? a : -a;
	}

	BinTree* findNearest(K k) {
		if(k == key) {
			return this;
		} else if(k < key) {
			if(!left) {
				return this;
			} else {
				BinTree* tmp = left->findNearest(k);
				if(abs(tmp->key - k) < abs(key - k)) {
					return tmp;
				} else {
					return this;
				}
			}
		} else if(k > key) {
			if(!right) {
				return this;
			} else {
				BinTree* tmp = right->findNearest(k);
				if(abs(tmp->key - k) < abs(key - k)) {
					return tmp;
				} else {
					return this;
				}
			}
		}
		return nullptr;
	}

public:
	K key;
	V value;

	BinTree(K k, V v) :
		key(k),
		value(v),
		left(nullptr),
		right(nullptr) {}

	BinTree* add(K k, V v) {
		if(k == key) {
			value = v;
			return this;
		} else if(k < key) {
			if(!left) {
				return (left = new BinTree(k, v));
			} else {
				return left->add(k, v);
			}
		} else if(k > key) {
			if(!right) {
				return (right = new BinTree(k, v));
			} else {
				return right->add(k, v);
			}
		}
		return nullptr;
	}

	bool findNearest(K k, K& actualK, V& v) {
		BinTree* n = findNearest(k);
		if(n) {
			actualK = n->key;
			v = n->value;
			return true;
		} else {
			return false;
		}
	}

	bool get(K k, V& v) {
		if(k == key) {
			v = value;
			return false;
		} else if(k < key) {
			return left ? left->get(k, v) : false;
		} else if(k > key) {
			return right ? right->get(k, v) : false;
		}
		return false;
	}

	~BinTree() {
		if(left) {
			delete left;
			left = nullptr;
		}
		if(right) {
			delete right;
			right = nullptr;
		}
	}
};

#endif /* INCLUDE_BINTREE_HPP_ */
