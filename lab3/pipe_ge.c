#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <mpi.h>
/*--------------- Gaussian elimination -- Pipelined version */
int pipe_ge(double *AA, int n,  MPI_Comm comm){
  /* Pipelined Gaussian Elimination  
     AA = pointer to matrix and right-hand side. AA is a one-dimensional 
	  array holding a matrix of size n*(n+1)-- row-wise. The
          right-hand side occupies the last column. 
     n  = dimension of problem
     On return Gaussian elimination has been applied and only upper triangular
     part is relevant -- back-solve should be called to get solution x.  
     *----------------------------------------------------------------------*/
  int nprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  int South = (myid+1)%nprocs;
  int North = (myid-1+nprocs)%nprocs;
  double* a_k = (double*)malloc(sizeof(double)*(n+1));
  double piv;
  /*-------------------- blank function + Do nothing and return */
  for (int k=0; k<n-1; k++){
    if (k%nprocs==myid){
      MPI_Send(&AA[(k/nprocs)*(n+1)], n+1, MPI_DOUBLE, South, k, comm);
      memcpy(a_k, &AA[(k/nprocs)*(n+1)], (n+1)*sizeof(double)); 
    }else{
      // If k is not in the last block, or k in the last block but myid is behind k%nprocs,
      // myid will need to use k.
      if ((k/nprocs<n/nprocs-1)||(k%nprocs<myid)) MPI_Recv(a_k, n+1, MPI_DOUBLE, North, k, comm, MPI_STATUS_IGNORE);
      // Same as above but check whether south has row k first.
      if (South!=(k % nprocs)&&((k/nprocs<n/nprocs-1)||(k%nprocs<South))) MPI_Send(a_k, n+1, MPI_DOUBLE, South, k, comm);
    }
    for (int i=k+1;i<n;i++){
      // check whether myid is in charge of row i
      if (i%nprocs==myid){
        piv = AA[i/nprocs*(n+1)+k]/a_k[k];
        for (int j=k+1;j<n+1;j++){
          AA[i/nprocs*(n+1)+j]-=piv*a_k[j];
        }
      }
    }
  }
  free(a_k);
  return(0);
}
