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
	edge_t *adj_list;
} vertex_t;

typedef struct {
	vertex_t *data;
	unsigned size;
} graph_t;

typedef struct {
	float distance;
	int vertex;
} weighted_vertex_t;

graph_t init_graph(unsigned size);
void delete_graph(graph_t g);

weighted_vertex_t* find_nearest_neighbours(graph_t g, int v, unsigned *listsz, float maxdist);

void link_edges(graph_t g, int v1, int v2, float dist);
void unlink_edge(graph_t g, int v);




#endif
