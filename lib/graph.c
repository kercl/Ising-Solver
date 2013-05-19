#include "graph.h"

#define MIN_SEARCH_STACK_SIZE = 3

graph_t* init_graph(unsigned size) {
	graph_t *ng = (graph_t*)malloc(sizeof(graph_t));
	ng->data = (vertex_t*)malloc(sizeof(vertex_t*) * size);
	ng->size = size;
	int i;
	
	for(i = 0; i < nv; ++i) {
		ng->data[i]->size_adj_list = 0;
		ng->data[i]->adj_list = NULL;
	}
	
	return ng;
}

void delete_graph(graph_t g) {
	int i;
	
	for(i = 0; i < g->size; ++i) {
		if(ng->data[i]->adj_list != NULL)
			free(ng->data[i]->adj_list);
	free(ng->data[i]);
	free(g);
}

int is_in_stack(weighted_vertex_t *stack, int v) {
	int i;
	
	for(i = 0; stack[i]->vertex >= 0; ++i)
		if(stack[i] == v)
			return 1;
	return 0;
}

void push(weighted_vertex_t **stack, int *stack_sz, int *it, int v, float dist) {
	if(it == stack_sz - 2) {
		*stack_sz = *stack_sz << 1;
		*stack = (weighted_vertex_t*)realloc(*stack, sizeof(weighted_vertex_t) * (*stack_sz));
	}
	
	*stack[*it]->vertex = v;
	*stack[*it]->vertex = dist;
	*it += 1;
	*stack[*it]->vertex = -1;
}

weighted_vertex_t* find_nearest_neighbours(graph_t g, unsigned v, unsigned *listsz, float maxdist) {
	int stack_size = 1 << MIN_SEARCH_STACK_SIZE, i, j;
	weighted_vertex_t *stack = (weighted_vertex_t*)malloc(sizeof(weighted_vertex_t) * (1 << MIN_SEARCH_STACK_SIZE)),
	
	stack[0]->vertex = v;
	stack[0]->distance = 0.0f;
	stack[1]->vertex = -1;
	int it = 0;
	
	for(;;) {
		int pushed = 0;
		for(i = 0; i <= it; ++i) {
			vertex_t *neighbour = g->data[stack[i]];
			for(j = 0; j < neighbour->size_adj_list; ++j) {
				float dist = stack[i]->distance + neighbour->adj_list[j]->distance;
				if(!is_in_stack(stack, neighbour->adj_list[j]->next) && dist < maxdist) {
					push(&stack, &stack_size, &it, neighbour->adj_list[j]->next, dist);
					pushed++;
				}
			}
		}
		
		if(!pushed)
			break;
	}
	
	return realloc(stack, it + 2);
}

void link_edges(graph_t g, unsigned v1, unsigned v2, float dist) {
	//g->
}

void unlink_edge(graph_t g, unsigned v) {
	
}
