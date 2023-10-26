CC = tcc 

TARGETS = *.c

LINKER_FLAGS = -g -o

OBJECT = ga

all: $(TARGETS)
	$(CC) $(TARGETS) $(LINKER_FLAGS) $(OBJECT) 
