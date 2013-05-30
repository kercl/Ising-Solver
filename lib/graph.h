#ifndef GRAPH_H
#define GRAPH_H

#include <stdint.h>
#include <math.h>

typedef struct {
	float distance;
	int a, b;
} edge_t;

typedef struct {
	unsigned size_adj_list;
	unsigned *adj;
	float *dist;
} vertex_t;

typedef struct {
	vertex_t *data;
	unsigned size;
	
	int *last_nn;  // result of the most recent nearest neighbour search
	float *distance; // the distance to the center of the sphere in the nn search
	unsigned last_nn_sz, last_nn_buffer_sz;
} graph_t;

typedef struct {
	float distance;
	int vertex;
} weighted_vertex_t;

void init_graph(graph_t *g, unsigned size);
void delete_graph(graph_t *g);

void graph_nearest_neighbours(graph_t *g, int v, float radius);

void link_edges(graph_t *g, unsigned v1, unsigned v2, float dist);
void unlink_edge(graph_t *g, unsigned v);

/* some graph generators */

void build_square_grid(graph_t *g, int w, int h);

#endif
