stencil: stencil.c
	gcc -std=c99 -pg -fopenmp -O3 -Wall $^ -o $@
