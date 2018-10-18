stencil: stencil.c
	gcc -std=c99 -O3 -pg -Wall $^ -o $@
