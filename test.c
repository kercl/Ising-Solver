#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lib/vector.h"
#include "lib/lattice.h"
#include "lib/ga.h"
#include "lib/graph.h"
#include "lib/avl.h"

#define DEFAULT_SPACING 15.0f

#define HIGHLIGHT 80
#define RADIUS 3.0f

#define WIDTH 8
#define HEIGHT 8

#define TYPE 0

void draw_dot(float x, float y, float r, int spin, FILE *fp) {
	char *col = "black";
	if(spin == 1)
		col = "red";
	else if(spin == -1)
		col = "blue";
	else if(spin == 2)
		col = "green";
	fprintf(fp, "<circle cx='%f' cy='%f' r='%f' stroke='%s' stroke-width='1' fill='%s'/>\n", x, y, r, col, col);
}

void draw_line(float x0, float y0, float x1, float y1, FILE *fp) {
	fprintf(fp, "<line x1='%f' y1='%f' x2='%f' y2='%f' style='stroke:rgb(0,0,0);stroke-width:0.05'/>", x0, y0, x1, y1);
}

float offset(int i, lattice_t *l) {
	switch(l->lattice_type) {
		case LATTICE_SQUARE_GRID:
			return 0;
		case LATTICE_TRIANGLE_GRID:
			return (i / l->width) & 1 ? 0 : DEFAULT_SPACING / 2.0f;
	}
	return 0.0f;
}

float offset_graph(int i, int type) {
	return (i / WIDTH) & 1 ? 0 : DEFAULT_SPACING / 2.0f;
}

void draw_grid(model_t *m, int psel) {
	FILE *fp = fopen("sim.svg", "w");
	
	fprintf(fp, "<svg xmlns='http://www.w3.org/2000/svg' version='1.1'>\n");
	
	if(m->topology_type == MODEL_TYPE_LATTICE) {
		float x0 = 1.0f, y0 = 1.0f;
		lattice_t *l = (lattice_t*)m->topology;
		for(int i = 0; i < l->width * l->height; ++i) {
			nearest_neighbours(l, i);
			
			if(l->last_nn_sz == 0)
				continue;
			
			draw_dot(x0 + offset(i, l) + (i % l->width) * DEFAULT_SPACING, 
					 y0 + (i / l->width) * DEFAULT_SPACING, 0.6, m->population_state[psel][i], fp);
			int j;
			for(j = 0; j < l->last_nn_sz; ++j)
				draw_line(x0 + offset(i, l) + (i % l->width) * DEFAULT_SPACING,
						  y0 + (i / l->width) * DEFAULT_SPACING,
						  x0 + offset(l->last_nn[j], l) + (l->last_nn[j] % l->width) * DEFAULT_SPACING,
						  y0 + (l->last_nn[j] / l->width) * DEFAULT_SPACING, fp);
		}
	}
	
	if(m->topology_type == MODEL_TYPE_GRAPH) {
		float x0 = 1.0f, y0 = 1.0f;
		graph_t *l = (graph_t*)m->topology;
		for(int i = 0; i < WIDTH * HEIGHT; ++i) {
			graph_nearest_neighbours(l, i, 1.42f);
			//getchar();
			
			if(l->last_nn_sz == 0)
				continue;
			
			draw_dot(x0 + offset_graph(i, TYPE) + (i % WIDTH) * DEFAULT_SPACING, 
					 y0 + (i / WIDTH) * DEFAULT_SPACING, 0.6, m->population_state[psel][i], fp);

			int j;
			for(j = 0; j < l->last_nn_sz; ++j)
				draw_line(x0 + offset_graph(i, TYPE) + (i % WIDTH) * DEFAULT_SPACING,
						  y0 + (i / WIDTH) * DEFAULT_SPACING,
						  x0 + offset_graph(l->last_nn[j], TYPE) + (l->last_nn[j] % WIDTH) * DEFAULT_SPACING,
						  y0 + (l->last_nn[j] / WIDTH) * DEFAULT_SPACING, fp);
		}
		int j;
		/*graph_nearest_neighbours(l, HIGHLIGHT, 1);

		int i;
		for(i = 0; i < l->last_nn_sz; ++i)
			draw_dot(x0 + offset_graph(l->last_nn[i], TYPE) + (l->last_nn[i] % WIDTH) * DEFAULT_SPACING, 
					 y0 + (l->last_nn[i] / WIDTH) * DEFAULT_SPACING, 1.0f, 2, fp);*/
	}
	fprintf(fp, "</svg>");
	fclose(fp);
}

int test_ga() {
	printf("Ising solver GA test\n");
	
	int state_set[] = {-1, 1}, i, j, k;
	int n;
	
	srand(time(NULL));
	
	model_t m;
	init_population(&m, MODEL_TYPE_GRAPH, 20, WIDTH*HEIGHT, state_set, 2);
	init_graph(m.topology, WIDTH * HEIGHT);
	
	build_triangle_grid(m.topology, WIDTH, HEIGHT);
	
	print_graph_connections(m.topology);
	
	draw_grid(&m, 0);
	
	printf("prebuilding edge list (press key to continue)\n");
	int t = clock();
	getchar();
	precalc_edge_list(&m);
	printf("finished: %d\n", (clock() - t));

	printf("SIZE CONNECTION LIST: %fMB\n", (m.connections * sizeof(edge_t)) / 1048576.0f);
	
	printf("simulating...\n");

	m.mutation_inhibitor = -0.2f;
	m.restrict_selection_percentage = 0.2f;

	/*printf("PREBUILD CONNECTION LIST:\n");
	for(j = 0; j < m.connections; ++j) {
		printf("%d ---> %d | %f\n", m.connection_list[j].a, m.connection_list[j].b, m.connection_list[j].distance);
	}*/
	getchar();
	for(j = 0; j < 1000; ++j) {
		//rate_fitness_ising(&m);
		/*for(i = 0; i < m.population_sz; ++i) {
			printf("%03d: ", i);
			for(k = 0; k < m.genomes; ++k)
				printf("%s", m.population_state[i][k] == -1 ? "-" : "+");
			printf(" | %f, %f\n", m.fitness[i], energy_ising(&m, i));
		}
		printf("\n");
		getchar();*/
		evolve(&m);
		//getchar();
	}

	
	return 0;
	//rate_fitness(&m);
	
	//draw_grid(&m, 0);
	printf("TEST");
	getchar();
	for(i = 0; i < 3; ++i);
	
	return 0;
}

void print_graph_connections(graph_t *g) {
	int i, j;
	for(i = 0; i < g->size; ++i) {
		for(j = 0; j < g->data[i].size_adj_list; ++j)
			printf("[%3d] ---> [%3d] (d: %f)\n", i, g->data[i].adj[j], g->data[i].dist[j]);
	}
}

int test_graph() {
	printf("Ising solver GA test\n");
	
	int state_set[] = {-1, 1}, i, j, k;
	int n;
	
	srand(time(NULL));
	
	model_t m;
	init_population(&m, MODEL_TYPE_GRAPH, 200, WIDTH*HEIGHT, state_set, 2);
	init_graph(m.topology, WIDTH * HEIGHT);
	
	printf("building graph");
	build_square_grid(m.topology, WIDTH, HEIGHT);
	
	print_graph_connections(m.topology);
	
	draw_grid(&m, 0);

	printf("pk\n");
	getchar();

	minimal_energies(&m, NULL, NULL);

	getchar();
	printf("GENETIC PART\n");
	
	m.mutation_inhibitor = -0.1f;
	m.restrict_selection_percentage = 0.2f;
	
	precalc_edge_list(&m);
	for(i = 0; i < 5000; ++i)
		evolve(&m);
	
	for(i = 0; i < m.population_sz; ++i) {
		printf("%03d: ", i);
		for(k = 0; k < m.genomes; ++k)
			printf("%s", m.population_state[i][k] == -1 ? "-" : "+");
		printf(" | %f, %f\n", m.fitness[i], energy_ising(&m, i));
	}
	
	return 0;
}

int main(int argc, char **argv) {
	return test_graph();
}
