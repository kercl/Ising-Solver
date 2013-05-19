#include <stdio.h>
#include <stdlib.h>

#include "lib/vector.h"
#include "lib/lattice.h"
#include "lib/ga.h"

#define DEFAULT_SPACING 15.0f

#define WIDTH 10
#define HEIGHT 10

void draw_dot(float x, float y, float r, int spin, FILE *fp) {
	char *col = "black";
	if(spin == 1)
		col = "red";
	else if(spin == -1)
		col = "blue";
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
	return 0;
}

void draw_grid(model_t *m, int psel) {
	FILE *fp = fopen("sim.svg", "w");
	
	fprintf(fp, "<svg xmlns='http://www.w3.org/2000/svg' version='1.1'>\n");
	
	if(m->topology_type == MODEL_TYPE_LATTICE) {
		float x0 = 1.0f, y0 = 1.0f;
		lattice_t *l = (lattice_t*)m->topology;
		for(int i = 0; i < l->width * l->height; ++i) {
			draw_dot(x0 + offset(i, l) + (i % l->width) * DEFAULT_SPACING, 
					 y0 + (i / l->width) * DEFAULT_SPACING, 0.6, m->population_state[psel][i], fp);
			nearest_neighbours(l, i);
			int j;
			for(j = 0; j < l->last_nn_sz; ++j)
				draw_line(x0 + offset(i, l) + (i % l->width) * DEFAULT_SPACING,
						  y0 + (i / l->width) * DEFAULT_SPACING,
						  x0 + offset(l->last_nn[j], l) + (l->last_nn[j] % l->width) * DEFAULT_SPACING,
						  y0 + (l->last_nn[j] / l->width) * DEFAULT_SPACING, fp);
		}
	}
	fprintf(fp, "</svg>");
	fclose(fp);
}

int main(int argc, char **argv) {
	printf("Ising solver test application\n");
	
	int state_set[] = {-1, 1}, i;
	
	srand(time(NULL));
	
	model_t m;
	init_population(&m, MODEL_TYPE_LATTICE, 10, WIDTH*HEIGHT, state_set, 2);
	init_lattice(m.topology, LATTICE_SQUARE_GRID, WIDTH, HEIGHT, 1);
	
	draw_grid(&m, 0);
	
	evolve(&m);
	for(i = 0; i < m.population_sz; ++i) {
		printf("Chromosome %d fitness: %f\n", i, m.fitness[i]);
	}
	
	return 0;
}
