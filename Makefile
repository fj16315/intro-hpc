stencilMPI: stencilMPI.c
	mpicc -std=c99 -pg -O3 -Wall $^ -o $@
