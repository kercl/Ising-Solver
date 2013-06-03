#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "lib/ga.h"
#include "lib/graph.h"
#include "lib/lattice.h"

#define WIDTH 8
#define HEIGHT 8

#define TEST_A_RUNS 50
#define TEST_A_GENERATIONS 5000

#define TEST_B_RUNS 50
#define TEST_B_GENERATIONS 5000

/***************************************
 * test selector
 **************************************/
#define PERFORM_B


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

void simulate(model_t *m, float mut_inh, float restr_brd, FILE *fp) {
#ifdef PERFORM_A
	#define TEST_RUNS TEST_A_RUNS
	#define TEST_GENERATIONS TEST_A_GENERATIONS
#endif

#ifdef PERFORM_B
	#define TEST_RUNS TEST_B_RUNS
	#define TEST_GENERATIONS TEST_B_GENERATIONS
#endif
	
	int i, j;
	
	m->mutation_inhibitor = mut_inh;
	m->restrict_selection_percentage = restr_brd;
	
	char output[TEST_A_RUNS][256];

	#pragma omp parallel for
	for(i = 0; i < TEST_RUNS; ++i) {
		randomise_population(m);
		for(j = 0; j < TEST_GENERATIONS; ++j) {
			evolve(m);

		}
		//print_population(m);
		//getchar();
		
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
	
	m.radius_of_influence = 1.0f;
	m.mutation_inhibitor = -0.2f;
	m.restrict_selection_percentage = 0.1f;
	
	const float mutation_inhs[] = {-0.6f, -0.4f, -0.2f, -0.0f, 0.2f, 0.4f, 0.6f},
				restrict_brds[] = {0.1f, 0.5f};
	
	FILE *fp = fopen("test_A/results", "w");
	
	fprintf(fp, "%d\t%d\t%d\t\%d\n\n", TEST_A_RUNS, TEST_A_GENERATIONS, WIDTH, HEIGHT);
	printf("Starting simulation test A\n");
	
	for(i = 0; i < sizeof(mutation_inhs) / sizeof(float); ++i)
		for(j = 0; j < sizeof(restrict_brds) / sizeof(float); ++j) {
			simulate(&m, mutation_inhs[i], restrict_brds[j], fp);
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
	 * DETERMINE ENERGY BY BRUTEFORCE
	 *************************************************************/
	 
#ifdef PERFORM_B

	init_graph(m.topology, WIDTH * HEIGHT);
	build_diamond_grid(m.topology, WIDTH, HEIGHT);
	
	m.radius_of_influence = 1.0f;
	
	printf("Starting simulation test B\n");

	FILE *fp = fopen("test_B/results", "w");

	simulate(&m, 0.6f, 0.3f, fp);
#endif
}
