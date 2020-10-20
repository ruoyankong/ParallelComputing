#include <stdlib.h> 
#include <stdio.h>
#include <mpi.h>
int back_solve(double *A, int n, int nloc, double *x, MPI_Comm comm){
  /*--------------------
  / upper triangular solution 
  / A = pointer to an n * (n+1) matrix. Only the upper triangular part of
        A is used. The right-hand side occupies the last column of A.
    The column version of the algorithm is implemented. 
  */
  /*-------------------- NOTE: THIS JUST SETS X = RHS AND RETURNS */
  /* it is provided so that you can see how to access cerain elements of
     the matrix stored as a one-dimensional array                 */
  int k, id,  nprocs, np=n+1, myid, kloc;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid); 
  double t;
  /*--------------- back-solve loop          */
  for (k=n-1; k>=0; k--) {
    id = (k % nprocs) ;
    if (id == myid) {
      kloc = (k-id)/nprocs;
/*-------------------- set x=rhs on return */ 
      t = A[kloc*np+n]/A[kloc*np+k];
      x[kloc] = t;
    }
    // barrier until all process updated and t does no change.
    MPI_Barrier(MPI_COMM_WORLD);
    // sent t to all the processes.
    MPI_Bcast(&t, 1, MPI_DOUBLE, id, comm);
    for (int i = 0; i<nloc; i++) A[i*np+n] -= t*A[i*np+k];
  }
  printf (" ------  end back-solve  in PE %d   ---------\n \n",myid);
  return(0);
}
