stencil: stencil.c
	icc -std=c99 -pg -O2 -Wall $^ -o $@
