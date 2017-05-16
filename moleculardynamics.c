#include <mpi.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <math.h>

#define FALSE 0
#define TRUE (!FALSE)

int abrt_on_err(int id, char* err_msg)
{
    if ( id == 0 ) {
        fprintf( stderr, err_msg );
    }
    MPI_Finalize();
    exit( 1 );
}

int balancer(int id, int p, int dflag)
{
    /*No balancing*/
    return p; 
    /*balancing*/
    if(dflag)
        return p;
    /*dflag is FALSE*/
    return 2 * (p - id) - 1;
}

partial_force_calc(double* X, int n, int p, int id, double* F)
{
    int i,j;
    int dflag;
    double tmp;

    for( i = 0; i < n; i++) {
        F[i] = 0.0;
    }

    dflag = TRUE;
    for( i = id; i < n - 1; i += balancer(id, p, dflag)) {
        for( j = i + 1; j < n; j++) {
            tmp = (double) 1.0 / ((X[i] - X[j]) * (X[i] - X[j]));
            F[i] += tmp;
            F[j] -= tmp;
        }
        dflag = !dflag;
    }
}


int main(int argc, char** argv) {
    int p, id, n;
    char *tp;
    int i, j;
    double *X, *F, *FF, tmp;

    /* Initialize the MPI environment */
    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Get the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if ( argc == 1 ) { /* No particle number given, getting out*/
        abrt_on_err(id, "No particle number given, aborting..\n");
    }
    errno = 0;
    n = strtol(argv[1], &tp, 10);
    if ( *tp != '\0' || errno != 0 ) { /*Checking for integer particle number*/
        abrt_on_err(id, "Wrong particle number, aborting..\n");
    }
    if (id == 0 ) { 
        printf("The number of particles is %ld\n", n);
    }
    if ( n > 20000 ) {
        abrt_on_err(id, "Particle number is too large (>20000), aborting..\n");
    }

    X = (double*)malloc(sizeof(double) * n);
    F = (double*)malloc(sizeof(double) * n);
    FF = (double*)malloc(sizeof(double) * n);

    if ( id == 0 ) { /*Can be parallelized with MPI_Allgather, sequental now*/
        srand(time(NULL));
	printf("Generating %ld random particle coordinates.\n", n);
        for( i = 0; i < n; i++) {
	    X[i] = rand() / (double)RAND_MAX; 
	    /* avoiding division by zero if two randoms generated equal */
            X[i] += (double)i * DBL_MIN;
	}
    }
    /*time step ierations must start here*/
    MPI_Bcast(X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    partial_force_calc(X, n, p, id, F);
    MPI_Allreduce(F, FF, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /*resulting forces are in FF array*/
    /*Particle motion computation ought to be here*/
    /*time step iteration must end here*/
    /* Finalize the MPI environment.*/
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
