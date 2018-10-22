stencil: stencil.c
	gcc -std=c99 -pg -O3 -Wall $^ -o $@
