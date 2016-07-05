#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// array holding our values
int a[1000000001];

int main(int argc, char *argv[]) {
	// start counting time
	double tstart = 0.0, ttaken;
	
	// check to make sure we have correct # of args
	if( argc != 3){
		printf("Incorrect number of arguments\n");
		exit(1);
	}

	// ceiling of prime nums to calc
	int n = atoi(argv[1]);
	// # of threads to use
	int t = atoi(argv[2]);

	// make sure N is > 2
	if(n <= 2) {
		printf("N must be > 2\n");
		exit(1);
	}
	// make sure t is a positive #
	if(t < 1) {
		printf("T must be a positive int\n");
		exit(1);
	}
	// set # of threads
	omp_set_num_threads(t);

	// start count for prime # gen algo
	tstart = omp_get_wtime();

	int i, j, counter, lastPrime;

	// iterate from 2 to ceiling of 1/2N
	for(i = 2; i <= (n+1)/2; i++) {
		if(!a[i]) { // if this # hasn't been marked as a multiple of another #
			a[i] = 1; 
			if(i <= sqrt(n)) { // if i exceeds sqrt(n), then we have marked all possible multiples as non-prime <= N
				#pragma omp parallel for 
				for(j = i * i; j <= n; j+=i) {
					a[j] = -1;
				}
			}
		}
	}
	// end count of prime # gen algo 
	ttaken = omp_get_wtime() - tstart;
	printf("Time take for the main part: %f\n", ttaken);

	// create file name
	char buffer[20];
	snprintf(buffer, 20, "%d%s", n, ".txt");

	// openm file to write to 
	FILE *f = fopen(buffer, "w");
	if(f == NULL) {
		printf("Error creating file!\n");
		exit(1);
	}

	counter = 1;
	lastPrime = 2;
	for(i = 2; i <= n; i++) {
		if(a[i] != -1) { // write to file
			fprintf(f, "%d, %d, %d\n", counter, i, (i - lastPrime));
			lastPrime = i;
			counter++;
		}
	}


	return 0;
}
