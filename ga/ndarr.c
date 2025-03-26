#include <stdlib.h>
#include <stdio.h>

void** ndarr(int* shape, int dims, size_t element_size)
{
	int i = 0;
	int k = 0;
	size_t total_elements = 1;
	size_t total_pointers = shape[dims-2];
	size_t total_length = 0;
	size_t total_size = 0;
	size_t level_size = 0;

	for(i = 0; i < dims; i++) {
		total_elements *= shape[i];
	}
	
	for(i = dims-3; i >= 0; i--) {
		total_pointers *= shape[i];
		total_pointers += shape[i];
	}
	
	total_size = total_elements * element_size + total_pointers * sizeof(void**);
	void** base_ptr = (void**) malloc(total_size);
	if (base_ptr == NULL)
	{
		exit(1);
	}

	void **base_data = base_ptr + total_pointers * sizeof(void**);
	
	void *data_iterator = base_data;
	i = total_pointers - 1;

	level_size = total_elements / shape[dims-1];
	k = total_pointers - level_size;
	int loop_counter = 0; /* debugging purposes, will be removed */
	/* set in-data pointers */
	while(k < total_pointers) {
		base_ptr[k] = data_iterator;
		loop_counter++;
		k++;
		data_iterator += element_size * shape[dims-1];
	}
	printf("data pointers assigned: %d\n", loop_counter);
	printf("data iterator increment: %d\n", element_size*shape[dims-1]);

	int start;
	int end = total_pointers;
	void *pointer_iterator;
	for(i = dims - 3; i >= 0; --i) {
		end -= level_size;
		level_size = level_size / shape[i+1];
		start = end - level_size; 
		pointer_iterator = base_ptr + sizeof(void*) * end;
		k = start;
		printf("start: %d\n", start);
		for(k = start; k < end; k++)
		{
			base_ptr[k] = pointer_iterator;
			pointer_iterator += shape[i+1] * sizeof(void*);
		}
	}
	printf("\n");
	printf("%ld\n", element_size);
	printf("%ld\n", sizeof(void*));
	return base_ptr;
}


int main(void)
{
	int shape[3] = {5, 4, 3};
	void*** ret = (void***) ndarr(shape, 3, sizeof(float));
	float* f = (float*) ret[0][0];
	f[2] = 3;
	return 0;
}
