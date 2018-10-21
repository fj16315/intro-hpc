stencil: stencil.c
	icc -std=c99 -O2 -pg -Wall $^ -o $@
