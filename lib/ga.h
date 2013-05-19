#ifndef GA_H
#define GA_H

#include "lattice.h"
#include "graph.h"

#define MODEL_TYPE_GRAPH 1
#define MODEL_TYPE_LATTICE 2

#define DEFAULT_POPULATION_SIZE 20

#define EVOLUTION_LIMIT 10

typedef struct {
	void *topology;
	int topology_type;
	
	int **population_state;
	float *fitness;
	
	int population_sz, individuals;
	
	edge_t *connection_list;
	unsigned connections;
}model_t;

void init_population(model_t *m, int tt, int population_size, int individuals, int *state_set, int states);
void prebuild_edgelist(model_t *m);

void evolve(model_t *m);

void rate_fitness_ising(model_t *m);

void selection(model_t *m);
void mix_genomes(model_t *m);
void mutate(model_t *m);

#endif
