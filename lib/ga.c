#include "ga.h"

#include <stdlib.h>

#define ABS(x) ( (x)>0 ? (x) : -(x) )
#define PARENT(x) ( (x)&1 ? ((x)-1)>>1 : ((x)-2)>>1 )

#define SELECTION_ELITIST_PERCENTAGE 0.1

void init_population(model_t *m, int tt, int population_size, int individuals, int *state_set, int states) {
	m->topology_type = tt;
	switch(tt) {
		case MODEL_TYPE_LATTICE:
			m->topology = malloc(sizeof(lattice_t));
			break;
	}
	m->population_sz = population_size;
	
	m->population_state = malloc(sizeof(int*) * population_size);
	m->fitness = malloc(sizeof(float) * population_size);
	
	m->individuals = individuals;
	m->connection_list = NULL;
	
	int i, j;
	for(i = 0; i < population_size; ++i) {
		m->population_state[i] = malloc(sizeof(int) * individuals);
		for(j = 0; j < individuals; ++j)
			m->population_state[i][j] = state_set[rand() % states];
	}
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
		for(i = 0; i < m->individuals; ++i) {
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
}

void rate_fitness_ising(model_t *m) {
	int i, j;
	
	float J = 1.0f, accum_fitness = 0.0f;
	
	for(i = 0; i < m->population_sz; ++i) {
		float H = 0;
		for(j = 0; j < m->connections; ++j)
			H += J * m->population_state[i][m->connection_list[j].a] * m->population_state[i][m->connection_list[j].b];
		m->fitness[i] = ABS(H);
		
		accum_fitness += m->fitness[i];
	}
	
	for(i = 0; i < m->population_sz; ++i)
		m->fitness[i] /= accum_fitness;
	
	sort_population(m);
}

void selection(model_t *m) {
	
}

void mix_genomes(model_t *m);
void mutate(model_t *m);
