stencil: stencil.c
	gcc -std=c99 -O2 -pg -Wall $^ -o $@
