#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>

#define BLOCK_SIZE 1024
#define WARP 32

long getmax(long *, long);
__global__ void get_block_max(long arr[], int size, int thread_size);

int main(int argc, char *argv[])
{
   long size = 0;  // The size of the array
   long i;  // loop index
   long * numbers; //pointer to the array
   long * array;
   long max = 0;
   
	cudaSetDevice(1);
    if(argc !=2)
    {
       printf("usage: maxseq num\n");
       printf("num = size of the array\n");
       exit(1);
    }

    size = atol(argv[1]);
    
	


    numbers = (long *)malloc(size * sizeof(long));
    if( !numbers )
    {
       printf("Unable to allocate mem for an array of size %ld\n", size);
       exit(1);
    }
	

	
    srand(time(NULL)); // setting a seed for the random number generator
    // Fill-up the array with random numbers from 0 to size-1
    for( i = 0; i < size; i++) {
       numbers[i] = rand() % size;
	}
	
	int elements = size - (size % WARP);
	int threads = (int)ceil((double)elements/BLOCK_SIZE);
	int blocks = (int)ceil(((double)threads/BLOCK_SIZE));
	
	cudaError_t err = cudaMalloc((void**)&array, sizeof(long) * elements);
	if (err != cudaSuccess) {
		printf("cudaMalloc failure\n");
	}
	err = cudaMemcpy(array, numbers, sizeof(long) * elements, cudaMemcpyHostToDevice);

	
	get_block_max<<<blocks, BLOCK_SIZE>>>(array, elements, BLOCK_SIZE);
	
	
	for (i = elements; i < size; i++) {
		if (max < numbers[i]) {
			max = numbers[i];
		}
	}
	
	cudaMemcpy(numbers, array, sizeof(long) * blocks, cudaMemcpyDeviceToHost);
	
	for (i = 0; i < blocks; i++) {
		if (numbers[i] > max) {
			max = numbers[i];
		}
	}

	
	
    printf(" The maximum number in the array is: %ld\n", max);

	cudaFree(array);
    free(numbers);
    exit(0);
}

// method that finds max within a block
__global__
void get_block_max(long arr[], int size, int thread_size) {
	int i;
	int max = 0;
	__shared__ int sdata[BLOCK_SIZE];
	
	unsigned int index = (threadIdx.x + blockIdx.x * blockDim.x) * thread_size;

	for (i = index; i < index + thread_size && i < size; i++) {
		if (max < arr[i]) {
			max = arr[i];
		}
	}
	sdata[threadIdx.x] = max;
	__syncthreads();

	unsigned int thid = threadIdx.x;
	unsigned int stride;

	for (stride = blockDim.x/2; stride >= WARP; stride >>= 1) {
		if (thid < stride && thid + stride < blockDim.x) {
			int curr = sdata[thid + stride];
			if (sdata[thid] < curr) {
				sdata[thid] = curr;
			}
		}
		__syncthreads();
	} 

	__syncthreads();

	if (thid == 0) {
		int i;
		int max = 0;

		for (i = thid; i < thid + WARP && i < size; i++) {
			if (max < sdata[i]) {
				max = sdata[i];
			}
		}
		sdata[0] = max;
		arr[blockIdx.x] = sdata[0];
	}
}




/*
   input: pointer to an array of long int
          number of elements in the array
   output: the maximum number of the array
*/
long getmax(long num[], long size)
{
  long i;
  long max = num[0];

  for(i = 1; i < size; i++)
	if(num[i] > max)
	   max = num[i];

  return( max );

}
