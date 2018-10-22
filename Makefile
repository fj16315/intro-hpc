stencil: stencil.c
	icc -std=c99 -pg -Ofast -Wall $^ -o $@
