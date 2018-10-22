stencil: stencil.c
	icc -std=c99 -pg -Wall $^ -o $@
