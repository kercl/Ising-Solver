#include "graph.h"
#include "avl.h"

#include <stdlib.h>

#define MIN_SEARCH_STACK_SIZE 7

#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

void init_graph(graph_t *ng, unsigned size) {
	ng->data = malloc(sizeof(vertex_t) * size);
	ng->size = size;
	int i;
	
	for(i = 0; i < size; ++i) {
		ng->data[i].size_adj_list = 0;
		ng->data[i].adj = NULL;
	}
	
	ng->last_nn = malloc((1 << MIN_SEARCH_STACK_SIZE) * sizeof(int));
	ng->distance = malloc((1 << MIN_SEARCH_STACK_SIZE) * sizeof(float));
	ng->last_nn_sz = 0;
	ng->last_nn_buffer_sz = 1 << MIN_SEARCH_STACK_SIZE;
}

void delete_graph(graph_t *g) {
	int i;
	
	for(i = 0; i < g->size; ++i) {
		if(g->data[i].adj != NULL)
			free(g->data[i].adj);
	}
	free(g->data);
	free(g);
}

void unwrap_tree(graph_t *g, avl_t *t) {
	if(!g || !t)
		return;
	
	g->last_nn[g->last_nn_sz] = t->a;
	g->distance[g->last_nn_sz] = t->dist;
	g->last_nn_sz++;
	
	unwrap_tree(g, t->l);
	unwrap_tree(g, t->r);
}

void graph_nearest_neighbours(graph_t *g, int v, float radius) {
	avl_t *searchtree = NULL;
	int newelement;
	
	g->last_nn_sz = 0;
	
	weighted_vertex_t *newbuffer = malloc((1 << MIN_SEARCH_STACK_SIZE) * sizeof(weighted_vertex_t)),
					  *searchbuffer = malloc((1 << MIN_SEARCH_STACK_SIZE) * sizeof(weighted_vertex_t));
	int size_search = 1, it_n = 0, i, j;
	
	searchbuffer[0].vertex = v;
	
	while(size_search > 0) {
		for(i = 0; i < size_search; ++i) {
			for(j = 0; j < g->data[searchbuffer[i].vertex].size_adj_list; ++j) {
				if(g->data[searchbuffer[i].vertex].adj[j] == v)
					continue;
				newelement = 0;
				
				if(g->data[searchbuffer[i].vertex].dist[j] + searchbuffer[i].distance <= radius)
					searchtree = insert(searchtree, g->data[searchbuffer[i].vertex].adj[j], 0,
													g->data[searchbuffer[i].vertex].dist[j] + searchbuffer[i].distance, 
													&newelement);
				if(newelement) {
					newbuffer[it_n].vertex = g->data[searchbuffer[i].vertex].adj[j];
					newbuffer[it_n].distance = g->data[searchbuffer[i].vertex].dist[j] + searchbuffer[i].distance;
					it_n++;
				}
			}
		}
		
		weighted_vertex_t *tmp = newbuffer;
		newbuffer = searchbuffer;
		searchbuffer = tmp;
		
		size_search = it_n;
		it_n = 0;
	}
	
	unwrap_tree(g, searchtree);
	delete_tree(searchtree);
}

void link_edge_unidir(graph_t *g, unsigned v1, unsigned v2, float dist) {
	int i;
	
	int updated = 0;
	for(i = 0; i < g->data[v1].size_adj_list; ++i)
		if(g->data[v1].adj[i] == v2) {
			g->data[v1].dist[i] = dist;
			updated = 1;
		}
	if(updated == 0) {
		g->data[v1].adj = realloc(g->data[v1].adj, (g->data[v1].size_adj_list + 1) * sizeof(unsigned));
		g->data[v1].dist = realloc(g->data[v1].dist, (g->data[v1].size_adj_list + 1) * sizeof(float));
		
		g->data[v1].adj[g->data[v1].size_adj_list] = v2;
		g->data[v1].dist[g->data[v1].size_adj_list] = dist;
		g->data[v1].size_adj_list++;
	}
}

void link_edges(graph_t *g, unsigned v1, unsigned v2, float dist) {	
	link_edge_unidir(g, v1, v2, dist);
	link_edge_unidir(g, v2, v1, dist);
}

void build_square_grid(graph_t *g, int w, int h) {
	int i;
	for(i = 0; i < w * h; ++i) {
		if(i % w + 1 < w) {
			link_edges(g, i, i+1, 1.0f);
		}
		if(i + w < w * h)
			link_edges(g, i, i+w, 1.0f);
	}
}

void build_triangle_grid(graph_t *g, int w, int h) {
	int i;
	for(i = 0; i < w * h; ++i) {
		if(i % w + 1 < w) {
			link_edges(g, i, i+1, 1.0f);
		}
		if(i + w < w * h) {
			link_edges(g, i, i+w, 1.0f);

			if((i / w) & 1) {
				if(i % w + 1 < w)
					link_edges(g, i, i+w+1, 1.0f);
			}else {
				if(i % w - 1 >= 0)
					link_edges(g, i, i+w-1, 1.0f);
			}
		}
	}
}

void build_diamond_grid(graph_t *g, int w, int h) {
	int i;
	
	for(i = 0; i < w * h; ++i) {
		if((i / w) & 1) { // odd row
			if( ((i % w) & 1) == 0 && (i % w) > 0) {
				link_edges(g, i, i-1, 1.0f);
			}
		}else { // even row
			if( ((i % w) & 1) == 0 ) { // x component
				if(i + w < w * h) {
					link_edges(g, i, i+w, 1.0f);
					if(i % w + 1 < w)
						link_edges(g, i, i+w+1, 1.0f);
				}
				if(i - w >= 0) {
					link_edges(g, i, i-w, 1.0f);
					if(i % w + 1 < w)
						link_edges(g, i, i-w+1, 1.0f);
				}
			} // otherwise, leave unconnected
		}
	}
}
