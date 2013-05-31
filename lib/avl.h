#ifndef AVL_H
#define AVL_H

struct avl {
	struct avl *l, *r;
	int a, b, dl, dr;
	float dist;
};

typedef struct avl avl_t;

avl_t* insert(avl_t *root, int a, int b, float dist, int *new);
void delete_tree(avl_t *t);

void debug_avl(avl_t *st, int intend, const char *prefix);

#endif
