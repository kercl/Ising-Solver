#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "lib/ga.h"
#include "lib/graph.h"
#include "lib/lattice.h"

#define DEFAULT_SPACING 15.0f

#define WIDTH 8
#define HEIGHT 8

#define TEST_A_RUNS 50
#define TEST_A_GENERATIONS 800

#define TEST_B_RUNS 30
#define TEST_B_GENERATIONS 800

#define TEST_C_RUNS 30
#define TEST_C_GENERATIONS 1000

/***************************************
 * test selector
 **************************************/
#define PERFORM_B

#ifdef PERFORM_C
	#define WIDTH 8
	#define HEIGHT 8
	#define TEST_C_GENERATIONS 1000
#endif

#ifdef PERFORM_D
	#define WIDTH 3
	#define HEIGHT 3
	#define TEST_C_GENERATIONS 1000
#endif

#define TEST_RUNS 0
#define TEST_GENERATIONS 0

#ifdef PERFORM_B
	int offset(int i, int subtest) {
		switch(subtest) {
			case 0:
				return 0;
			case 1:
				return (i / WIDTH) & 1 ? DEFAULT_SPACING / 2.0f : 0;
			case 2:
				return (i / WIDTH) & 1 ? 0 : DEFAULT_SPACING / 2.0f;
		}
		return 0;
	}
#endif

#ifdef PERFORM_A
	int offset(int i, int subtest) {
		return 0;
	}
#endif

#ifdef PERFORM_D
	int offset(int i, int subtest) {
		return 0;
	}
#endif

void draw_dot(float x, float y, float r, int spin, FILE *fp) {
	char *col = "black";
	if(spin == 1)
		col = "black";
	else if(spin == -1)
		col = "white";
	else if(spin == 2)
		col = "green";
	fprintf(fp, "<circle cx='%f' cy='%f' r='%f' stroke='black' stroke-width='0.3' fill='%s'/>\n", x, y, r, col, col);
}

void draw_line(float x0, float y0, float x1, float y1, float sw, FILE *fp) {
	fprintf(fp, "<line x1='%f' y1='%f' x2='%f' y2='%f' style='stroke:rgb(0,0,0);stroke-width:%f'/>", x0, y0, x1, y1, sw);
}

void print_graph_connections(graph_t *g) {
	int i, j;
	for(i = 0; i < g->size; ++i) {
		printf("NODE: %d - ", i);
		for(j = 0; j < g->data[i].size_adj_list; ++j)
			printf("%d[%f]  ", g->data[i].adj[j], g->data[i].dist[j]);
		printf("\n");
	}	
}

void draw_graph(model_t *m, char *fn, int psel, int subtest) {
	FILE *fp = fopen(fn, "w");
	
	fprintf(fp, "<svg xmlns='http://www.w3.org/2000/svg' version='1.1'>\n");
	
	if(m->topology_type == MODEL_TYPE_GRAPH) {
		float x0 = 1.0f, y0 = 1.0f;
		graph_t *l = (graph_t*)m->topology;
		for(int i = 0; i < WIDTH * HEIGHT; ++i) {
			graph_nearest_neighbours(l, i, 1.42f);
			
			if(l->last_nn_sz == 0)
				continue;
			
			draw_dot(x0 + offset(i, subtest) + (i % WIDTH) * DEFAULT_SPACING, 
					 y0 + (i / WIDTH) * DEFAULT_SPACING, 1.2, m->population_state[psel][i], fp);

			int j;
			for(j = 0; j < l->last_nn_sz; ++j)
				draw_line(x0 + offset(i, subtest) + (i % WIDTH) * DEFAULT_SPACING,
						  y0 + (i / WIDTH) * DEFAULT_SPACING,
						  x0 + offset(l->last_nn[j], subtest) + (l->last_nn[j] % WIDTH) * DEFAULT_SPACING,
						  y0 + (l->last_nn[j] / WIDTH) * DEFAULT_SPACING, 
						  0.1 / (2 * l->distance[j]),
						  fp);
		}
	}
	fprintf(fp, "</svg>");
	fclose(fp);
}

int spinsum(model_t *m, int i) {
	int s = 0, j;
	for(j = 0; j < m->genomes; ++j) {
		s += m->population_state[i][j];
	}
	return s;
}

void print_population(model_t *m) {
	int i, k;
	for(i = 0; i < m->population_sz; ++i) {
		printf("%03d: \n", i);
		for(k = 0; k < m->genomes; ++k) {
			if(k>0 && k%100==0)
				printf("\n");
			printf("%s", m->population_state[i][k] == -1 ? "-" : "+");
		}
		printf(" | %f, %f\n", m->fitness[i], energy_ising(m, i));
	}
	printf("\n");
}

void simulate(model_t *m, float mut_inh, float restr_brd, FILE *fp, int subtest) {
#ifdef PERFORM_A
	#define TEST_RUNS TEST_A_RUNS
	#define TEST_GENERATIONS TEST_A_GENERATIONS
	
	char pic_prefix[100];
	sprintf(pic_prefix, "test_A/sim_convergence_%06d_", subtest);
#endif

#ifdef PERFORM_B
	#define TEST_RUNS TEST_B_RUNS
	#define TEST_GENERATIONS TEST_B_GENERATIONS
	
	char pic_prefix[100];
	sprintf(pic_prefix, "test_B/sim_pic_%d_", subtest);
#endif
	
	int i, j;
	
	m->mutation_inhibitor = mut_inh;
	m->restrict_selection_percentage = restr_brd;
	
	char output[TEST_A_RUNS][256];

	#pragma omp parallel for
	for(i = 0; i < TEST_RUNS; ++i) {
		randomise_population(m);
	#ifdef PERFORM_A
		char tmp[100];
		sprintf(tmp, "%s_%d_graph", pic_prefix, i); 
		FILE *fp_converge = fopen(tmp, "w");
	#endif
		
		for(j = 0; j < TEST_GENERATIONS; ++j) {
		#ifdef PERFORM_A
			rate_fitness(m);
			fprintf(fp_converge, "%d\t%f\n", j, energy_ising(m, 0));
		#endif
			evolve(m);
		}
	#ifdef PERFORM_A
		fclose(fp_converge);
	#endif
		//print_population(m);
		//draw_graph(m, "tmptmp.svg", 0, 0);
		//getchar();
		
		rate_fitness(m);

	#ifdef PERFORM_B
		char tmp[100];
		sprintf(tmp, "%s_%d.svg", pic_prefix, i);
		draw_graph(m, tmp, 0, subtest);
	#endif
		
		printf("\rRun %3d/%3d completed  ", i+1, TEST_RUNS);
		fflush(stdout);
		
		sprintf(output[i], "%f\t%f\t%f\t%d\n", mut_inh, restr_brd, energy_ising(m, 0), spinsum(m, 0));
	}
	for(i = 0; i < TEST_RUNS; ++i) {
		fprintf(fp, output[i]);
	}
	fflush(fp);
}

int main(int argc, char **argv) {
	int state_set[] = {-1, 1};
	int i, j, k;
	
	char str[100];
	
	model_t m;
	init_population(&m, MODEL_TYPE_GRAPH, 200, WIDTH*HEIGHT, state_set, 2);
	
	
	/*************************************************************
	 * TEST A:
	 * TESTING OF THE GENETIC ALGORITHM WITH THE BASIC ISING MODEL
	 *************************************************************/
	 
#ifdef PERFORM_A

	init_graph(m.topology, WIDTH * HEIGHT);
	build_square_grid(m.topology, WIDTH, HEIGHT);
	
	m.radius_of_influence = 1.1f;
	m.mutation_inhibitor = -0.2f;
	m.restrict_selection_percentage = 0.1f;
	
	const float mutation_inhs[] = {-5.0f, -4.8f, -4.6f, -4.4f, -4.2f},
				restrict_brds[] = {0.1f};
	
	FILE *fp = fopen("test_A/results", "w");
	
	fprintf(fp, "%d\t%d\t%d\t\%d\n\n", TEST_A_RUNS, TEST_A_GENERATIONS, WIDTH, HEIGHT);
	printf("Starting simulation test A\n");
	
	precalc_edge_list(&m);
	printf("%d\n", m.connections);
	for(i = 0; i < sizeof(mutation_inhs) / sizeof(float); ++i)
		for(j = 0; j < sizeof(restrict_brds) / sizeof(float); ++j) {
			simulate(&m, mutation_inhs[i], restrict_brds[j], fp, i*10000+j);
			printf("\nFinished configuration %d/%d\n", i*(sizeof(restrict_brds) / sizeof(float))+j+1, (sizeof(mutation_inhs) / sizeof(float))*(sizeof(restrict_brds) / sizeof(float)));
		}
	
	fclose(fp);

#endif

	/*************************************************************
	 * TEST B:
	 * GEOMETRIC FRUSTRATION ON LATTICES
	 *  _    .     |
	 * [_]  /_\  _/ \ _
	 *            \ /
	 *             |
	 *************************************************************/
	 
#ifdef PERFORM_B
	printf("Starting simulation test B\n");
	m.radius_of_influence = 1.1f;

	init_graph(m.topology, WIDTH * HEIGHT);
	build_square_grid(m.topology, WIDTH, HEIGHT);
	
	precalc_edge_list(&m);
	m.J = malloc(m.connections * sizeof(float));
	for(i = 0; i < m.connections; ++i)
		m.J[i] = -1.0f;
	
	FILE *fp = fopen("test_B/results_0", "w");

	simulate(&m, -3.2f, 0.1f, fp, 0);
	
	fclose(fp);
	delete_graph(m.topology);
	
	init_graph(m.topology, WIDTH * HEIGHT);
	build_triangle_grid(m.topology, WIDTH, HEIGHT);

	precalc_edge_list(&m);
	free(m.J);
	m.J = malloc(m.connections * sizeof(float));
	for(i = 0; i < m.connections; ++i)
		m.J[i] = -1.0f;
	
	fp = fopen("test_B/results_1", "w");

	simulate(&m, -3.2f, 0.1f, fp, 1);

	fclose(fp);
	delete_graph(m.topology);
	
	init_graph(m.topology, WIDTH * HEIGHT);
	build_diamond_grid(m.topology, WIDTH, HEIGHT);
	
	precalc_edge_list(&m);
	free(m.J);
	m.J = malloc(m.connections * sizeof(float));
	for(i = 0; i < m.connections; ++i)
		m.J[i] = -1.0f;
	
	fp = fopen("test_B/results_2", "w");

	simulate(&m, -3.2f, 0.1f, fp, 2);
	
	fclose(fp);
#endif

	/*************************************************************
	 * TEST C:
	 * SPIN GLASS ON TRIANGULAR 8x8 GRID
	 *************************************************************/
	 
#ifdef PERFORM_C
	delete_model(&m);
	init_population(&m, MODEL_TYPE_GRAPH, 300, WIDTH*HEIGHT, state_set, 2);

	init_graph(m.topology, WIDTH * HEIGHT);
	build_triangle_grid(m.topology, WIDTH, HEIGHT);
	
	precalc_edge_list(&m);
	m.J = malloc(m.connections * sizeof(float));
	for(i = 0; i < m.connections; ++i)
		m.J[i] = random01() * 2.0f - 1.0f;
	
	m.radius_of_influence = 1.1f;
	
	FILE *fp = fopen("test_C/results", "w");
	
	fprintf(fp, "%d\t%d\t%d\t\%d\n\n", TEST_C_RUNS, TEST_C_GENERATIONS, WIDTH, HEIGHT);
	printf("Starting simulation test C\n");
	
	simulate(&m, 0.6f, 0.05f, fp, 0);
	
	fclose(fp);

#endif

#ifdef PERFORM_D
	delete_model(&m);
	init_population(&m, MODEL_TYPE_GRAPH, 300, WIDTH*HEIGHT, state_set, 2);

	init_graph(m.topology, WIDTH * HEIGHT);
	build_triangle_grid(m.topology, WIDTH, HEIGHT);
	
	precalc_edge_list(&m);
	
	minimal_energies(&m, NULL, NULL);
#endif
}
