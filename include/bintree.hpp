/*
 * bintree.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: rob
 */

#ifndef INCLUDE_BINTREE_HPP_
#define INCLUDE_BINTREE_HPP_

namespace hlrg {
namespace ds {

/**
 * This is a simple binary tree which can be treated like a map, or
 * used for nearest-match searches.
 */
template <class K, class V>
class BinTree {
private:
	BinTree<K, V>* left;
	BinTree<K, V>* right;

	/**
	 * Return the absolut value of a.
	 *
	 * @param a A numeric value.
	 * @return The absolute value of the input.
	 */
	inline K abs(K a) {
		return a > 0 ? a : -a;
	}

	/**
	 * Return the node containing the nearest match to k.
	 *
	 * @param k A numeric key.
	 * @return  The node nearest the key, or a null pointer if nothing was found, though that should be impossible.
	 */
	BinTree* findNearest(K k) {
		if(k == key) {
			// If the key is equal to this node's key, return this node.
			return this;
		} else if(k < key) {
			if(!left) {
				// If there's no left branch, return this node.
				return this;
			} else {
				// Find the nearest branch on the left to k.
				BinTree* tmp = left->findNearest(k);
				if(tmp) {
					if(abs(tmp->key - k) < abs(key - k)) {
						// If it is closer than this node's key, return it.
						return tmp;
					} else {
						// Otherwise return this node.
						return this;
					}
				}
			}
		} else if(k > key) {
			if(!right) {
				// If there's no right branch, return this node.
				return this;
			} else {
				// Find the nearest branch on the right to k.
				BinTree* tmp = right->findNearest(k);
				if(tmp) {
					if(abs(tmp->key - k) < abs(key - k)) {
						// If it is closer than this node's key, return it.
						return tmp;
					} else {
						// Otherwise return this node.
						return this;
					}
				}
			}
		}
		// It should be a failure to arrive here.
		return nullptr;
	}

public:
	K key;
	V value;

	/**
	 * Construct a BinTree using the given key and value.
	 *
	 * @param k A key value.
	 * @param v A value.
	 */
	BinTree(K k, V v) :
		left(nullptr),
		right(nullptr),
		key(k),
		value(v) {}

	/**
	 * Add the given key and value to the tree. If the key already exists,
	 * the value will be overwritten.
	 *
	 * @param k A key value.
	 * @param v A value.
	 */
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

	/**
	 * Find the nearest key and value to the given key. The actualK argument
	 * will be updated with the value of the nearest key, and the v argument
	 * will receive the value.
	 *
	 * Returns true on success, false if nothing is found.
	 *
	 * @param k The search key.
	 * @param actualK A reference to a key which will be updated with the closest key.
	 * @param v A reference to a value that will be updated with the closest value.
	 * @return True if the search succeeds.
	 */
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

	/**
	 * Find the exact match for the given key and update the value.
	 *
	 * Returns true if the value is found, false otherwise.
	 *
	 * @param k A search key.
	 * @param v A value to update.
	 * @return True if the key is found, false otherwise.
	 */
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

} // ds
} // hlrg

#endif /* INCLUDE_BINTREE_HPP_ */
