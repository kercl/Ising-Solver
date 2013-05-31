#include "avl.h"

#include <stdlib.h>

#define AVL_R 2
#define AVL_L 1
#define ROTATE_L 1
#define ROTATE_R 2

#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

void debug_avl(avl_t *st, int intend, const char *prefix) {

	if(!st) {
		return;
	}
	int i = 0;
	for(i = 0; i < intend; ++i)
		printf("  ");
	printf("%s[%d, %d, %f | %d, %d] <= %X\n", prefix, st->a, st->b, st->dist, st->dl, st->dr, (int)st & 0xFFFF);
	debug_avl(st->l, intend+1, "l");
	debug_avl(st->r, intend+1, "r");
}

/*****************************************************
 * Free memory, occupied by the tree
 *****************************************************/

void delete_tree(avl_t *t) {
	if(!t)
		return;
	delete_tree(t->l);
	delete_tree(t->r);
	
	free(t);
}

/*****************************************************
 * Perform a rotation of a subtree
 * adapts the depth values of the tree automatically
 * dir ... direction: either ROTATE_L or ROTATE_R
 *****************************************************/

avl_t* rotated(avl_t *root, int dir) {
	if(dir == ROTATE_R) {
		avl_t *ret = root->l,
			  *nd = ret->r;
		ret->r = root;
		root->l = nd;
		
		ret->r->dl = ret->dr;
		ret->dr = MAX(ret->r->dl, ret->r->dr) + 1;
		return ret;
	}else {
		avl_t *ret = root->r,
			  *nd = ret->l;
		ret->l = root;
		root->r = nd;
		
		ret->l->dr = ret->dl;
		ret->dl = MAX(ret->l->dr, ret->l->dl) + 1;
		return ret;
	}
}

/*****************************************************
 * Balances a subtree root according to the rules
 * of AVL trees
 *****************************************************/

avl_t* balance(avl_t *root) {
	if(root->dl - root->dr == 2) {
		if(root->l->dr > root->l->dl)
			root->l = rotated(root->l, ROTATE_L);
		return rotated(root, ROTATE_R);
	}
	else if(root->dr - root->dl == 2) {
		if(root->r->dl > root->r->dr)
			root->r = rotated(root->r, ROTATE_R);
		return rotated(root, ROTATE_L);
	}
	return root;
}

/*****************************************************
 * Insert operation 
 * root: the tree where the new element sould be
 * 		 added to.
 * a: index a of the new element
 * b: index b of the new element
 * dist: float, user data
 * new: after completion, the variable referenced by
 * 		new contains ether 1, if a new element was
 * 		added or 0 otherwise
 * returns the updated tree
 *****************************************************/

avl_t* insert(avl_t *root, int a, int b, float dist, int *new) {
	int n;
	if(new == NULL)
		new = &n;
	
	// if the tree is still empty, create a new element.
	if(root == NULL) {
		avl_t *n = malloc(sizeof(avl_t));
		n->a = a;
		n->b = b;
		n->dist = dist;
		n->l = NULL;
		n->r = NULL;
		n->dl = 0;
		n->dr = 0;
		*new = 1;
		return n;
	}
	
	*new = 0;
	
	// stack used to propagate back to the root element, without having to use recursive function
	avl_t **stack = malloc((MAX(root->dl, root->dr) + 2)* sizeof(avl_t*));
	int spos = 1;
	stack[0] = root;
	
	avl_t *it = root, *rot;
	int nextl = 0;
	for(;;) {
		nextl = 0;
		stack[spos++] = it;
		
		// decide, whether the new element should be inserted in the left or right subtree
		if(it->a < a)
			nextl = AVL_R;
		else if(it->a > a)
			nextl = AVL_L;
		else {
			if(it->b < b)
				nextl = AVL_R;
			else if(it->b > b)
				nextl = AVL_L;
		}
		
		if(nextl == 0) {
			free(stack);
			return root; // node already exists
		}
		
		// if the next layer is empty, create new element and abort, otherwise sink to next layer and repead the loop
		if(nextl == AVL_R) {
			if(it->r != NULL)
				it = it->r;
			else {
				it->r = insert(NULL, a, b, dist, new);
				it = it->r;
				break;
			}
		}else if(nextl == AVL_L) {
			if(it->l != NULL)
				it = it->l;
			else {
				it->l = insert(NULL, a, b, dist, new);
				it = it->l;
				break;
			}
		}
	}
	
	// adjust depth information and balance tree
	*new = 1;
	for(spos--; spos >= 0; --spos) {
		int dl = 0, dr = 0;
		
		if(it == stack[spos]->l && stack[spos]->l != NULL) {
			stack[spos]->l = balance(stack[spos]->l);
			stack[spos]->dl = MAX(stack[spos]->l->dl, stack[spos]->l->dr) + 1;
		}
		else if(it == stack[spos]->r && stack[spos]->r != NULL) {
			stack[spos]->r = balance(stack[spos]->r);
			stack[spos]->dr = MAX(stack[spos]->r->dl, stack[spos]->r->dr) + 1;
		}
		it = stack[spos];
	}
	
	free(stack);
	
	return balance(root);
}
