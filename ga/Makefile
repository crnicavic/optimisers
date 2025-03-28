CC = gcc

LINKER_FLAGS = -g -lm -o -Wall


lib: genetic.c genetic.h
	cc -fPIC -shared -o libgenetic.so genetic.c

test: test.c genetic.c genetic.h
	$(CC) genetic.c test.c $(LINKER_FLAGS) test
