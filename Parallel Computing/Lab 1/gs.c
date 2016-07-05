#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/***** Globals ******/
float *a; /* The coefficients */
float *x;  /* The unknowns */
float *b;  /* The constants */
float *curr; /* current values of unknown */
float err; /* The absolute relative error */
int num;  /* number of unknowns */
int MPI_Init(int *argc, char ***argv);

/****** Function declarations */
void check_matrix(); /* Check whether the matrix will converge */
void get_input();  /* Read input from file */
int calc_error();

/********************************/



/* Function definitions: functions are ordered alphabetically ****/
/*****************************************************************/

/*
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
  int bigger = 0; /* Set to 1 if at least one diag element > sum  */
  int i, j;
  float sum = 0;
  float aii = 0;


  for(i = 0; i < num; i++)
  {
    sum = 0;
    aii = fabs(a[i * num + i]);

//    printf("%d\n", a[i][i]);

    for(j = 0; j < num; j++)
       if( j != i)
	 sum += fabs(a[i * num + j]);

    if( aii < sum)
    {
      printf("The matrix will not converge\n");
      exit(1);
    }

    if(aii > sum)
      bigger++;

  }

  if( !bigger )
  {
     printf("The matrix will not converge\n");
     exit(1);
  }
}


/******************************************************/
/* Read input from file */
void get_input(char filename[])
{
  FILE * fp;
  int i,j;

  fp = fopen(filename, "r");
  if(!fp)
  {
    printf("Cannot open file %s\n", filename);
    exit(1);
  }

 fscanf(fp,"%d ",&num);
 fscanf(fp,"%f ",&err);

// printf("%d\n", num);

 /* Now, time to allocate the matrices and vectors */
 a = (float*)malloc(num * num * sizeof(float));
 if( !a)
  {
	printf("Cannot allocate a!\n");
	exit(1);
  }

 x = (float *) malloc(num * sizeof(float));
 if( !x)
  {
	printf("Cannot allocate x!\n");
	exit(1);
  }


 b = (float *) malloc(num * sizeof(float));
 if( !b)
  {
	printf("Cannot allocate b!\n");
	exit(1);
  }

 /* Now .. Filling the blanks */

 /* The initial values of Xs */
 for(i = 0; i < num; i++)
	fscanf(fp,"%f ", &x[i]);


 for(i = 0; i < num; i++)
 {
   for(j = 0; j < num; j++)
     fscanf(fp,"%f ",&a[i * num + j]);
   /* reading the b element */
 fscanf(fp,"%f ",&b[i]);
 }

 fclose(fp);

}

int calc_error(float* cur_x, int num) {
	int i;
	float error;
	int flag = 0;
	for(i = 0; i < num; i++) {
		error = fabs((cur_x[i] - x[i]) / cur_x[i]);
//		printf("%f\n", error);
		if (error > err){
			flag = 1;
			i = num;
		}
	}
	return flag;
}


/************************************************************/


int main(int argc, char *argv[])
{

 int i, j, k;
 int rank, size;
 int nit = 0; /* number of iterations */



 if( argc != 2)
 {
   printf("Usage: gsref filename\n");
   exit(1);
 }

 /* Read the input file and fill the global data structure above */
 get_input(argv[1]);

 /* Check for convergence condition */
 check_matrix();


 // start up MPI
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);





 float* cur_x = (float *)malloc(num * sizeof(float)); // current x for calc
 int* sizes = (int *)malloc(size * sizeof(int)); // # of rows each process receives
 int* sizes_rows = (int *)malloc(size * sizeof(int)); // # of elements each process receives for a

 int* start_index = (int *)malloc(size * sizeof(int)); // starting row index
 int* start_index_rows = (int *)malloc(size * sizeof(int)); // starting index for elem in a
 start_index[0] = 0; // initalize to 0 for process 0

 int rows = num / size; // min # of rows
 int extra_rows = num % size;// # of extra rows to allocate to earlier processes

 for(i = 0; i < size; i++) {
	 sizes[i] = rows;
	 if (i < extra_rows) sizes[i]++;// add an extra row if the process is low enough #
	 sizes_rows[i] = sizes[i] * num;

	 if (i > 0) start_index[i] = start_index[i-1] + sizes[i-1];
	 start_index_rows[i] = start_index[i] * num;
 }


 // initialize local vars for each process
 //float* loc_a = (float *)malloc(sizes_rows[rank] * sizeof(float));
 float* loc_x = (float *)malloc(sizes[rank] * sizeof(float));
 //float* loc_b = (float *)malloc(sizes[rank] * sizeof(float));

//scatter a, b, and x for calculations
 //MPI_Scatterv(a, sizes_rows, start_index_rows, MPI_FLOAT, loc_a, sizes_rows[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);
 //MPI_Scatterv(b, sizes, start_index, MPI_FLOAT, loc_b, sizes[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);
// MPI_Scatterv(x, sizes, start_index, MPI_FLOAT, loc_x, sizes[rank], MPI_FLOAT, 0, MPI_COMM_WORLD);

// initalize cur x
 for(i = 0; i < num; i++){
     cur_x[i] = x[i];
 }
do{
	nit++;
    for(i = 0; i < num; i++){
        x[i] = cur_x[i];
    }
	for(i = 0; i < sizes[rank]; i++) {
		 float sum = 0;
		 for(j = 0; j < num; j++) {
			 if(start_index[rank] + i == j){ // skip diag
				 continue;
			 }
			 else{
				 sum = sum + a[start_index_rows[rank] + i * num + j] * x[j]; // summate x[n] * a[n]
			 }
		 }
		 loc_x[i] = (b[start_index[rank] + i] - sum) / a[start_index_rows[rank] + i * num + start_index[rank] + i]; // set local x using jacobi equation
	}
	MPI_Allgatherv(loc_x, sizes[rank], MPI_FLOAT, cur_x, sizes, start_index, MPI_FLOAT, MPI_COMM_WORLD); // gather cur x from each process and update cur x
}while(calc_error(cur_x, num));



 /* Writing to the stdout */
 /* Keep that same format */
if(rank == 0) {
	for( i = 0; i < num; i++)
	printf("%f\n",cur_x[i]);
	printf("total number of iterations: %d\n", nit);
}



// MPI_Barrier(MPI_COMM_WORLD);
 MPI_Finalize();


 exit(0);

}
