#include <stdlib.h>
#include <stdio.h>

#include "lib/ga.h"
#include "lib/graph.h"
#include "lib/lattice.h"

#define WIDTH 10
#define HEIGHT 10

int spinsum(model_t *m, int i) {
	int s = 0, j;
	for(j = 0; j < m->genomes; ++j) {
		s += m->population_state[i][j];
	}
	return s;
}

int main(int argc, char **argv) {
	int state_set[] = {-1, 1};
	int i, j, k;
	
	char str[100];
	
	model_t m;
	init_population(&m, MODEL_TYPE_GRAPH, 20, WIDTH*HEIGHT, state_set, 2);
	
	init_graph(m.topology, WIDTH * HEIGHT);
	build_square_grid(m.topology, WIDTH, HEIGHT);
	
	m.radius_of_influence = 1.0f;
	m.mutation_inhibitor = -0.2f;
	m.restrict_selection_percentage = 0.1f;
	
	/*************************************************************
	 * TEST A:
	 * TESTING OF THE GENETIC ALGORITHM WITH THE BASIC ISING MODEL
	 * CONVERGENCE BEHAVIOUR
	 *************************************************************/
	 
	#define TEST_A_RUNS 5
	#define TEST_A_GENERATIONS 10000
	
	FILE *fp_m = fopen("test_A_data", "w");
	
	fprintf(fp_m, "runs: %d\ngenerations: %d\n", TEST_A_RUNS, TEST_A_GENERATIONS);
	printf("Starting simulation test A\n");
	for(i = 0; i < TEST_A_RUNS; ++i) {
		sprintf(str, "test_A_convergence.frm%02d", i);
		FILE *fp = fopen(str, "w");

		randomise_population(&m);
		for(j = 0; j < TEST_A_GENERATIONS; ++j) {
			evolve(&m);
			rate_fitness(&m);
			fprintf(fp, "%04d %f\n", j, energy_ising(&m, 0));
		}
		
		fprintf(fp_m, "> gen %d:\n  best case energy: %f\n  spin sum: %d\n\n", i, energy_ising(&m, 0), spinsum(&m, 0));
		
		fclose(fp);
	}
	fclose(fp_m);
	
	/*************************************************************
	 * TEST B:
	 * GEOMETRIC FRUSTRATION ON LATTICES
	 *  _    .     |
	 * [_], /_\, _/ \ _
	 *            \ /
	 *             |
	 * DETERMINE ENERGY BY BRUTEFORCE
	 * Edwards–Anderson model
	 *************************************************************/
	 
	/*************************************************************
	 * 
}
