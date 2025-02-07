CC = gcc

EX1 = genetic.c example1.c

EX2 = genetic.c example2.c

LINKER_FLAGS = -g -lm -o

OBJECT_EX1 = ex1

OBJECT_EX2 = ex2

lib: genetic.c genetic.h
	cc -fPIC -shared -o libgenetic.so genetic.c

example1: example1.c genetic.c genetic.h
	$(CC) $(EX1) $(LINKER_FLAGS) $(OBJECT_EX1) 

example2: example2.c genetic.c genetic.h
	$(CC) $(EX2) $(LINKER_FLAGS) $(OBJECT_EX2)

