#include <stdlib.h>
#include <stdio.h>

void* ndarr(int* shape, int dims, size_t element_size)
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

	void **end_ptr = base_ptr + total_size;
	float **base_data = (float**) base_ptr + total_pointers * sizeof(void**);
	
	void **data_iterator = (void**) base_data;
	i = total_pointers - 1;

	level_size = total_elements / shape[dims-1];
	/* set in-data pointers */
	while(k >= (total_pointers - level_size) && data_iterator < end_ptr) {
		base_ptr[k] = data_iterator;
		k--;
		data_iterator += element_size * shape[dims-1];
	}

	int start;
	int end;
	for(i = dims - 2; i > 0; --i) {
		level_size = level_size / shape[i];
		start = level_size / shape[(i-1)]; 
		end = start + level_size;
		printf("%d\n", end - start);
	}

	return NULL;
}


int main(void)
{
	int shape[4] = {6, 5, 4, 3};
	void *p1 = shape;
	void *p2 = &p1;
	void *p3 = &p2;
	int ***p3_int = (int***) p3;
	int a = p3_int[0][0][1];
	void* ret = ndarr(shape, 4, sizeof(float));
	return 0;
}
