CC = gcc 

TARGETS = *.h *.c

LINKER_FLAGS = -lm -g -o

OBJECT = ga

all: $(TARGETS)
	$(CC) $(TARGETS) $(LINKER_FLAGS) $(OBJECT)
	rm genetic.h.gch
