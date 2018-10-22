stencil: stencil.c
	icc -std=c99 -pg -O3 -Wall $^ -o $@
