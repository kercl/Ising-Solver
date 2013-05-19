#ifndef LATTICE_H
#define LATTICE_H

#define LATTICE_ARB 0
#define LATTICE_SQUARE_GRID 1
#define LATTICE_TRIANGLE_GRID 2

#include "vector.h"

typedef struct {
	int *nearest_displacements, nearest_displacements_sz;
	int (* disturb_fct)(int i, int displ, int w, int h, int d);
	
	int width, height, depth;
	
	int lattice_type;
	
	int *defect_list; /* only manipulate with add_defects and remove_defects */
	float *weights;
	
	unsigned defect_list_sz, weights_sz;
	
	int *last_nn;
	unsigned last_nn_sz, last_nn_buffer_sz;
}lattice_t;

void init_lattice(lattice_t *m, int type, int width, int height, int depth);

void add_defects(lattice_t *df, int *defects, int len);
void remove_defects(lattice_t *df, int *defects, int len);

void set_weights(lattice_t *df, i_vector3_t **edges, float *weights, int len);

float get_weight(lattice_t *df, i_vector3_t *a, i_vector3_t *b);

void nearest_neighbours_weight(lattice_t *df, int p, float rad);
void nearest_neighbours(lattice_t *df, int p);

#endif
