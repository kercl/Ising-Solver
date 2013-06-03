#include "ga.h"
#include "avl.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#define M_PI 3.14159265358979323846

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) c: (b) )

#define ABS(x) ( (x)>0 ? (x) : -(x) )
#define PARENT(x) ( (x)&1 ? ((x)-1)>>1 : ((x)-2)>>1 )

float random01() {
	return ((float)(rand() & 0xFFFF) / 65536.0f);
}

double gaussrand() {
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);

	phase = 1 - phase;

	return Z;
}

/*****************************************************
 * Initialise new population
 * m: model data structure holding the informations
 * tt: type of the topology
 * 		can be either MODEL_TYPE_LATTICE or
 * 		MODEL_TYPE_GRAPH
 * population_size: the number of individuals
 * genomes: the length of the characteristic 
 * 			chromosome string
 * state_set: possible allele states (typically {0,1})
 * states: number of possible alleles
 *****************************************************/
 
void init_population(model_t *m, int tt, int population_size, int genomes, int *state_set, int states) {
	m->topology_type = tt;
	switch(tt) {
		case MODEL_TYPE_LATTICE:
			m->topology = malloc(sizeof(lattice_t));
			break;
		case MODEL_TYPE_GRAPH:
			m->topology = malloc(sizeof(graph_t));
			break;
	}
	m->population_sz = population_size;
	
	m->population_state = malloc(sizeof(int*) * population_size);
	m->fitness = malloc(sizeof(float) * population_size);
	
	m->genomes = genomes;
	m->connection_list = NULL;
	m->connections = 0;
	
	m->state_set = state_set;
	m->n_states = states;
	
	m->restrict_selection_percentage = 1.0f;
	m->radius_of_influence = 1.0f;
	
	int i;
	for(i = 0; i < population_size; ++i)
		m->population_state[i] = malloc(sizeof(int) * m->genomes);
	
	randomise_population(m);
}

/*****************************************************
 * Initialise the genomes with random values
 *****************************************************/

void randomise_population(model_t *m) {
	int i, j;
	for(i = 0; i < m->population_sz; ++i) {
		for(j = 0; j < m->genomes; ++j) {
			m->population_state[i][j] = m->state_set[rand() % m->n_states];
		}
	}
}

/*****************************************************
 * Free occupied space
 *****************************************************/
 
void delete_model(model_t *m) {
	free(m->connection_list);
	
	int i;
	for(i = 0; i < m->population_sz; ++i)
		free(m->population_state[i]);
	free(m->population_state);
	
	if(m->topology_type == MODEL_TYPE_GRAPH)
		delete_graph(m->topology);
}

/*****************************************************
 * Completely delete a population
 *****************************************************/
 
void die_out(model_t *m) {
	int i;
	for(i = 0; i < m->population_sz; ++i)
		free(m->population_state[i]);
	free(m->population_state);
}

/*****************************************************
 * Expand the content of an AVL tree, containing
 * all possible connection information in the 
 * graph or lattice and store it in the connection
 * prebuild list (much faster computation time,
 * but very, very greedy considering memory
 * usage
 *****************************************************/
 
void unwrap_connections(model_t *g, avl_t *t) {
	if(!g || !t)
		return;
	
	g->connection_list[g->connections].a = t->a;
	g->connection_list[g->connections].b = t->b;
	g->connection_list[g->connections].distance = t->dist;
	g->connections++;
	
	unwrap_connections(g, t->l);
	unwrap_connections(g, t->r);
}

/*****************************************************
 * Prebuild a list of all possible pairs in a graph
 * where an interaction can occur
 * algorithm build an AVL tree (insertation O(log(n)))
 * instead of O(n), since a connection could be added
 * twice otherwise (from a to b and b to a)
 * the tree content can then easily be converted to
 * an array (unwrap_connections)
 *****************************************************/
 
void precalc_edge_list(model_t *m) {
	int i, j, k, ea, eb;
		
	if(m->connection_list != NULL) {
		free(m->connection_list);
		m->connection_list = NULL;
		m->connections = 0;
	}
	
	int edges = 0, new;
	avl_t *stree = NULL;
	
	if(m->topology_type == MODEL_TYPE_LATTICE) {
		lattice_t *l = m->topology;
		
		for(i = 0; i < m->genomes; ++i) {
			nearest_neighbours(l, i);
			
			for(k = 0; k < l->last_nn_sz; ++k) {
				if(i < l->last_nn[k]) {
					ea = i;
					eb = l->last_nn[k];
				}else {
					ea = l->last_nn[k];
					eb = i;
				}
				
				stree = insert(stree, ea, eb, l->distance[k], &new);
				if(new)
					edges++;
			}
		}
		
		m->connections = 0;
		m->connection_list = malloc(edges * sizeof(edge_t));
		unwrap_connections(m, stree);
		delete_tree(stree);
	}
	
	if(m->topology_type == MODEL_TYPE_GRAPH) {
		
		graph_t *l = m->topology;
		
		if(m->connection_list != NULL) {
			free(m->connection_list);
			m->connection_list = NULL;
		}
			
		for(i = 0; i < m->genomes; ++i) {
			graph_nearest_neighbours(l, i, m->radius_of_influence);
			
			for(k = 0; k < l->last_nn_sz; ++k) {
				if(i < l->last_nn[k]) {
					ea = i;
					eb = l->last_nn[k];
				}else {
					ea = l->last_nn[k];
					eb = i;
				}
				
				stree = insert(stree, ea, eb, l->distance[k], &new);
				if(new)
					edges++;
			}
		}
		
		m->connections = 0;
		m->connection_list = malloc(edges * sizeof(edge_t));
		unwrap_connections(m, stree);
		delete_tree(stree);
	}
}

/*****************************************************
 * Exchanges to populations
 *****************************************************/
 
void swap_population(model_t *m, int a, int b) {
	float f = m->fitness[a];
	m->fitness[a] = m->fitness[b];
	m->fitness[b] = f;
	
	int *p = m->population_state[a];
	m->population_state[a] = m->population_state[b];
	m->population_state[b] = p;
}

/*****************************************************
 * Rearrange populations, such that the list of
 * populations can be interpreted as a heap
 *****************************************************/
 
void heapify(model_t *m) {
	int i;
	
	float *f = m->fitness;
	
	for(i = 1; i < m->population_sz; ++i) {
		int pullback = i, parent = PARENT(i);
		while(f[parent] > f[pullback]) {
			swap_population(m, parent, pullback);
			pullback = parent;
			parent = PARENT(pullback);
		}
	}
}

/*****************************************************
 * Sort population according to their fitness
 * Algorithm: Heapsort
 *****************************************************/
 
void sort_population(model_t *m) {
	heapify(m);
	
	int end = m->population_sz - 1, pos = 0;
	do {
		swap_population(m, end, 0);
		end--;
		
		pos = 0;
		while((pos<<1) + 1 <= end) {
			if(m->fitness[(pos<<1)+1] < m->fitness[pos]) {
				if(m->fitness[(pos<<1)+2] < m->fitness[(pos<<1)+1] && (pos << 1)+2 <= end) {
					swap_population(m, (pos<<1)+2, pos);
					pos = (pos<<1)+2;
				}else {
					swap_population(m, (pos<<1)+1, pos);
					pos = (pos<<1)+1;
				}
			}else if(m->fitness[(pos<<1)+2] < m->fitness[pos] && (pos << 1)+2 <= end) {
				swap_population(m, (pos<<1)+2, pos);
				pos = (pos<<1)+2;
			}else
				break;
		}
	}while(end > 1);
	swap_population(m, 0, 1);
}

/*****************************************************
 * Breed a new generation
 *****************************************************/
 
void evolve(model_t *m) {
	if(m->connection_list == NULL)
		precalc_edge_list(m);
	if(m->connection_list == NULL)
		return;
	
	rate_fitness(m);
	breed(m);
	mutate(m);
}

/*****************************************************
 * Simple energy in the Ising model:
 *     ___
 *    \
 * H = )   J   S  S
 *    /___  ij  i  j
 *    <i,j>
 * 
 * J   = J constant
 *  ij
 *****************************************************/
 
float energy_ising(struct model *m, int i) {
	float H = 0.0f, J = 1.0f;
	int j;
	
	for(j = 0; j < m->connections; ++j)
		H += J * m->population_state[i][m->connection_list[j].a] * m->population_state[i][m->connection_list[j].b];
	return H;
}

/*****************************************************
 * Rates the fitness of each individual in a population
 *****************************************************/
 
void rate_fitness(model_t *m) {
	int i, j;
	
	float accum_fitness = 0.0f, minimum = INFINITY;
	
	for(i = 0; i < m->population_sz; ++i) {
		float H = energy_ising(m, i);
		m->fitness[i] = H;
		
		if(H < minimum)
			minimum = H;
	}
	
	for(i = 0; i < m->population_sz; ++i) {
		m->fitness[i] += ABS(minimum);
		accum_fitness += m->fitness[i];
	}
	
	for(i = 0; i < m->population_sz; ++i)
		m->fitness[i] /= accum_fitness;
	
	sort_population(m);
}

/*****************************************************
 * Randomly select an individual
 *****************************************************/

int select_random_individual(model_t *m, float accumfit) {
	float fit = random01() * accumfit;
	int i;
	
	float fitsum = 0.0f;
	for(i = 0; i < m->population_sz; ++i) {
		fitsum += m->fitness[i];
		if(fit <= fitsum)
			return i;
	}
	return 0;
}

/*****************************************************
 * Create a new generation (crossover)
 *****************************************************/

void breed(model_t *m) {
	int n = (int)(m->population_sz * m->restrict_selection_percentage), i, j, k;
	float accum_fitness = 0.0f;
	if(n == 0)
		n = 1;
	if(n >= m->population_sz) {
		n = m->population_sz;
		accum_fitness = 1.0f;
	}else {
		for(i = 0; i < n; ++i)
			accum_fitness += m->fitness[i];
	}
	
	/*
	 * Selection phase
	 */
	
	int **new_population = malloc(m->population_sz * sizeof(int*));
	
	for(i = 0; i < m->population_sz; ++i) {
		int indiv_a = select_random_individual(m, accum_fitness), 
			indiv_b = select_random_individual(m, accum_fitness);
		if(indiv_b == indiv_a) { // cannot select the same individual for breeding, prefer fitter individuals
			if(indiv_a == 0)
				indiv_b = 1;
			else
				indiv_b--;
		}
		
		int splice_at = rand() % m->genomes;

		// crossover at splice_at
		new_population[i] = malloc(m->genomes * sizeof(int));
		memcpy(new_population[i], m->population_state[indiv_a], splice_at * sizeof(int));
		memcpy(new_population[i] + splice_at, m->population_state[indiv_b] + splice_at, (m->genomes - splice_at) * sizeof(int));
	}
	
	// replace old population
	die_out(m);
	m->population_state = new_population;
}

/*****************************************************
 * Randomly flip bits in the individual's chromosomes
 *****************************************************/
 
void mutate(model_t *m) {
	int i, j;
	
	for(i = 0; i < m->population_sz; ++i)
		for(j = 0; j < m->genomes; ++j) {
			if(random01() < (1.0f - m->mutation_inhibitor) / ((float)(m->genomes)) )
				m->population_state[i][j] = m->state_set[rand() % m->n_states];
		}
}

/*****************************************************
 * Find minimal enery configurations by bruteforce
 *****************************************************/
 
void nth_state_to_model(model_t *m, uint8_t *nstate) {  /* currently only 2-state models !!! */
	int allele = 0;
	
	while(allele < m->genomes) {
		m->population_state[0][allele] = m->state_set[nstate[allele]];
		allele++;
	}
}

void inc_state(model_t *m, uint8_t *state) {
	int i = 0, carry = 1;
	while(carry) {
		state[i]++;
		carry = state[i] >= m->n_states;
		state[i] = state[i] % m->n_states;
		i++;
	}
}

int finished(model_t *m) {
	int i;
	for(i = 0; i < m->genomes; ++i)
		if(m->population_state[0][i] != m->state_set[m->n_states - 1])
			return 0;
	return 1;
}

void minimal_energies(model_t *template, int ***min_configs, int *n) {
	model_t m;
	init_population(&m, template->topology_type, 1, template->genomes, template->state_set, template->n_states);
	
	m.topology = template->topology;
	
	precalc_edge_list(&m);
	
	float running_max = INFINITY;
	
	uint8_t *state = malloc(template->genomes + 1);
	
	memset(state, 0, template->genomes * template->n_states);
	
	int k;
	while(!finished(&m)) {
		nth_state_to_model(&m, state);
		/*for(k = 0; k < m.genomes; ++k)
			printf("%s", m.population_state[0][k] == -1 ? "-" : "+");
		printf(" | %f\r", energy_ising(&m, 0));*/
		inc_state(&m, state);
	}
}
