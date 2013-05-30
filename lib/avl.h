#ifndef AVL_H
#define AVL_H

struct avl {
	struct avl *l, *r;
	int a;
	float dist;
};

typedef struct avl avl_t;

avl_t* insert(avl_t *root, int a, float dist, int *new);

#endif
