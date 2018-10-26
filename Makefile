stencil: stencil.c
	icc -std=c99 -pg -vec -xHost -Ofast -Wall $^ -o $@
