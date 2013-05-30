#include "avl.h"

struct avl {
	struct avl *l, *r;
	int a;
	float dist;
};

typedef struct avl avl_t;

avl_t* rotated(avl_t *root, int dir) {
	if(dir == ROTATE_R) {
		avl_t *ret = root->l,
			  *nd = ret->r;
		ret->r = root;
		root->l = nd;
		return ret;
	}else {
		avl_t *ret = root->r,
			  *nd = ret->l;
		ret->l = root;
		root->r = nd;
		return ret;
	}
}

avl_t* rebalance(avl_t *root, int *h) {
	int hl, hr;
	
	if(root == NULL) {
		*h = 0;
		return NULL;
	}
	
	root->l = rebalance(root->l, &hl);
	root->r = rebalance(root->r, &hr);
	
	if(hl - hr == 2) {
		*h = (hl + hr) / 2 + 1;
		return rotated(root, ROTATE_R);
	}
	if(hr - hl == 2) {
		*h = (hl + hr) / 2 + 1;
		return rotated(root, ROTATE_L);
	}
	
	*h = MAX(hl, hr) + 1;
	
	return root;
}

avl_t* insert(avl_t *root, int a, float dist, int *new) {
	if(root == NULL) {
		avl_t *n = malloc(sizeof(avl_t));
		n->a = a;
		n->dist = dist;
		n->l = NULL;
		n->r = NULL;
		*new = 1;
		return n;
	}
	
	*new = 0;
	
	avl_t *it = root;
	int nextl = 0;
	for(;;) {
		nextl = 0;
		
		if(it->a < a)
			nextl = AVL_R;
		if(it->a > a)
			nextl = AVL_L;
		
		if(nextl == 0) 
			return root; // node already exists
		
		if(nextl == AVL_R) {
			if(it->r != NULL)
				it = it->r;
			else {
				it->r = insert(NULL, a, dist, new);
				break;
			}
		}
		if(nextl == AVL_L) {
			if(it->l != NULL)
				it = it->l;
			else {
				it->l = insert(NULL, a, dist, new);
				break;
			}
		}
	}
	
	*new = 1;
	
	int h;
	return rebalance(root, &h);
}
