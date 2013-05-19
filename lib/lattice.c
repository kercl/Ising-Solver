#include "lattice.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define CANTOR(x,y) ( (((x+y)*(x+y+1))>>1)+y )

void sort(int *data, int len) {
	int i, j, min, tmp;
	
	for(i = 0; i < len; ++i) {
		min = i;
		for(j = i + 1; j < len; ++j) {
			if(data[j] < data[tmp])
				tmp = j;
		}
		if(min != i) {
			tmp = data[min];
			data[min] = data[i];
			data[i] = tmp;
		}
	}
}

int disturbance_triangle(int i, int displ, int w, int h, int d) {
	if(displ <= -w)
		return 1-((i / h) & 1);
	if(displ >= w)
		return -((i / h) & 1);
	return 0;
}

void init_lattice(lattice_t *m, int type, int width, int height, int depth) {
	m->last_nn = malloc(32 * sizeof(unsigned));
	m->last_nn_buffer_sz = 32;
	m->last_nn_sz = 0;
	
	m->width = width;
	m->height = height;
	m->depth = depth;
	
	m->disturb_fct = NULL;
	
	switch(type) {
		case LATTICE_SQUARE_GRID:
			m->nearest_displacements_sz = 4;
			m->nearest_displacements = malloc(4 * sizeof(int));
			m->nearest_displacements[0] = 1;
			m->nearest_displacements[1] = -width;
			m->nearest_displacements[2] = -1;
			m->nearest_displacements[3] = width;
			break;
		case LATTICE_TRIANGLE_GRID:
			m->nearest_displacements_sz = 6;
			m->nearest_displacements = malloc(5 * sizeof(int));
			m->nearest_displacements[0] = 1;
			m->nearest_displacements[1] = -width;
			m->nearest_displacements[2] = -width-1;
			m->nearest_displacements[3] = -1;
			m->nearest_displacements[4] = width;
			m->nearest_displacements[5] = width+1;
			m->disturb_fct = &disturbance_triangle;
			break;
	}
	
	m->lattice_type = type;
}

void add_defects(lattice_t *df, int *defect, int len) {
	df->defect_list = realloc(df->defect_list, (df->defect_list_sz + len)*sizeof(int));
	memcpy(df->defect_list + df->defect_list_sz, defect, sizeof(int)*len);
	df->defect_list_sz += len;
	sort(df->defect_list, df->defect_list_sz);
}

int is_defect(lattice_t *df, int i) {  /* binary search for i in the defects list */
	int splitlen = df->defect_list_sz >> 1, splitpos = splitlen;
	while(splitlen > 0) {
		if(df->defect_list[splitpos] == i)
			return 1;
		
		splitlen >>= 1;
		if(df->defect_list[splitpos] > i)
			splitpos -= splitlen;
		else
			splitpos += splitlen;
	}
	
	return 0;
}

void remove_defects(lattice_t *df, int *defects, int len) {
	/* TODO */
}

void set_weights(lattice_t *df, i_vector3_t **edges, float *weights, int len) {
	
}

float get_weight(lattice_t *df, i_vector3_t *a, i_vector3_t *b) {	
	int pos = CANTOR( CANTOR(a->x,CANTOR(a->y,a->z)), CANTOR(b->x,CANTOR(b->y,b->z)) );
	
	if(pos < df->weights_sz)
		return df->weights[pos] == NAN ? 1.0f : df->weights[pos];
	
	return 1.0f;
}

void nearest_neighbours_weight(lattice_t *df, int p, float rad) {
	
}

void nearest_neighbours(lattice_t *df, int p) {
	df->last_nn_sz = 0;
	int i, gridsize = df->width * df->depth * df->height, j = 0;
	
	int dist_border_left = p % df->width,
		dist_border_right = p % df->width - (df->width - 1);
	
	for(i = 0; i < df->nearest_displacements_sz; ++i) {
		int neighbour = p + df->nearest_displacements[i];
		
		if(df->disturb_fct)
			neighbour += df->disturb_fct(p, df->nearest_displacements[i], df->width, df->height, df->depth);
		
		if(dist_border_left + df->nearest_displacements[i] % df->width < 0 ||
			dist_border_right + df->nearest_displacements[i] % df->width > 0)
			continue;
		
		if(neighbour >= 0 && neighbour < gridsize) {
			df->last_nn[j] = neighbour;
			j++;
		}
	}
	df->last_nn_sz = j;
}
