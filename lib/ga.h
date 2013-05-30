#ifndef GA_H
#define GA_H

#include "lattice.h"
#include "graph.h"

#define MODEL_TYPE_GRAPH 1
#define MODEL_TYPE_LATTICE 2

#define DEFAULT_POPULATION_SIZE 20

#define EVOLUTION_LIMIT 10
#define SELECTION_ELITIST_PERCENTAGE 0.2

typedef struct {
	void *topology;
	int topology_type;
	
	int **population_state;
	float *fitness;
	
	int population_sz, genomes;
	
	float radius_of_influence;
	
	edge_t *connection_list;
	unsigned connections;
	
	float restrict_selection_percentage;
	
	int *state_set;
	int n_states;
	
	float mutation_inhibitor;
}model_t;

void init_population(model_t *m, int tt, int population_size, int genomes, int *state_set, int states);
void precalc_edge_list(model_t *m);

void evolve(model_t *m);

void rate_fitness_ising(model_t *m);
float energy_ising(model_t *m, int i);

void breed(model_t *m);
void mutate(model_t *m);

#endif
