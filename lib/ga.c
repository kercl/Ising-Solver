#include "ga.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

#define ABS(x) ( (x)>0 ? (x) : -(x) )
#define PARENT(x) ( (x)&1 ? ((x)-1)>>1 : ((x)-2)>>1 )

float random01() {
	return ((float)(rand() & 0xFFFF) / 65536.0f);
}

void init_population(model_t *m, int tt, int population_size, int genomes, int *state_set, int states) {
	m->topology_type = tt;
	switch(tt) {
		case MODEL_TYPE_LATTICE:
			m->topology = malloc(sizeof(lattice_t));
			break;
	}
	m->population_sz = population_size;
	
	m->population_state = malloc(sizeof(int*) * population_size);
	m->fitness = malloc(sizeof(float) * population_size);
	
	m->genomes = genomes;
	m->connection_list = NULL;
	
	int i, j;
	for(i = 0; i < population_size; ++i) {
		m->population_state[i] = malloc(sizeof(int) * genomes);
		for(j = 0; j < genomes; ++j)
			m->population_state[i][j] = state_set[rand() % states];
	}
	
	m->restrict_selection_percentage = 1.0f;
	m->state_set = state_set;
	m->n_states = states;
}

void die_out(model_t *m) {
	int i;
	for(i = 0; i < m->population_sz; ++i)
		free(m->population_state[i]);
	free(m->population_state);
}

void precalc_edge_list(model_t *m) {
	int i, j, k, ea, eb;
		
	if(m->connection_list != NULL) {
		free(m->connection_list);
		m->connection_list = malloc(sizeof(edge_t));
	}
	
	int edges = 0;
	
	lattice_t *l = m->topology;
	
	if(m->topology_type == MODEL_TYPE_LATTICE) {
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
				
				for(j = 0; j < edges; ++j) {
					if(m->connection_list[j].a == ea && m->connection_list[j].b == eb)
						break;
				}
				if(j >= edges) {
					m->connection_list = realloc(m->connection_list, (edges + 1) * sizeof(edge_t));
					m->connection_list[edges].a = ea;
					m->connection_list[edges].b = eb;
					m->connection_list[edges].distance = 1.0f;
					edges++;
				}
			}
		}
	}
	
	m->connections = edges;
}

void swap_population(model_t *m, int a, int b) {
	float f = m->fitness[a];
	m->fitness[a] = m->fitness[b];
	m->fitness[b] = f;
	
	int *p = m->population_state[a];
	m->population_state[a] = m->population_state[b];
	m->population_state[b] = p;
}

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

void evolve(model_t *m) {
	if(m->connection_list == NULL)
		precalc_edge_list(m);
	
	rate_fitness_ising(m);
	breed(m);
	mutate(m);
}

float energy_ising(model_t *m, int i) {
	float H = 0.0f, J = -1.0f;
	int j;
	
	for(j = 0; j < m->connections; ++j)
		H += J * m->population_state[i][m->connection_list[j].a] * m->population_state[i][m->connection_list[j].b];
	return H;
}

void rate_fitness_ising(model_t *m) {
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

int* clone_population(model_t *m, int p) {
	int *tmp = malloc(m->genomes * sizeof(int));
	memcpy(tmp, m->population_state[p], m->genomes * sizeof(int));
	return tmp;
}

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
		while(indiv_b == indiv_a)
			indiv_b = select_random_individual(m, accum_fitness);
		
		int splice_at = rand() % m->genomes;

		new_population[i] = malloc(m->genomes * sizeof(int));
		memcpy(new_population[i], m->population_state[indiv_a], splice_at * sizeof(int));
		memcpy(new_population[i] + splice_at, m->population_state[indiv_b] + splice_at, (m->genomes - splice_at) * sizeof(int));
	}
	
	die_out(m);
	m->population_state = new_population;
}

void mutate(model_t *m) {
	int i, j;
	
	for(i = 0; i < m->population_sz; ++i)
		for(j = 0; j < m->genomes; ++j) {
			if(random01() < (1.0f - m->mutation_inhibitor) / ((float)(m->genomes)) )
				m->population_state[i][j] = m->state_set[rand() % m->n_states];
		}
}
