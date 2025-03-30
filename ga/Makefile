CC = gcc

LINKER_FLAGS = -lm -Wall -o

DEBUG_FLAGS = -g -lm -O0 -Wall -o

libgenetic.so: genetic.c genetic.h
	cc -fPIC -shared -o libgenetic.so genetic.c

test: test.c genetic.c genetic.h
	$(CC) genetic.c test.c $(LINKER_FLAGS) test

debug: test.c genetic.c genetic.h
	$(CC) genetic.c test.c $(LINKER_FLAGS) debug
	
clean:
	rm libgenetic.so
	rm test
	rm debug
