#include <stdlib.h> 
#include <string.h> 
#include <stdio.h> 
#include <mpi.h>
#include <math.h>

int HQuicksort(int *Abuf, int *lenA,  int *list, int N, MPI_Comm comm) {
/*--------------------    
 Most of  the lab2  assignment will consist  of writing  this function.
 Right now,  it is set  to just print ``nothing  done" and exit.  It is
 provided so you can compile and run the pgm - but the numbers will not
 be sorted.
-------------------- */
  int n, tmp;
  n = *lenA;
  int size = N / n;
  if (ceil(log2(size)) != floor(log2(size))) return(1);
  int myid;
  int ncub = log2(size);
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  int j, q, i, l, r, piv, ncount;
  int *number_buf, *b;
  for (j=ncub; j>0; j--){
    q = pow(2, j-1);
    i = pow(2, ncub-j)-1+myid/(2*q);
    piv = list[i];

    l = 0;
    r = n;
    while(l<r){
      while(Abuf[l] <= piv && l < r) l++;
      while(Abuf[r] > piv && l < r) r--;
      if (l >= r) break;
      tmp = Abuf[l];
      Abuf[l] = Abuf[r];
      Abuf[r] = tmp;
    }
    if (Abuf[l] > piv) l--;

    if (myid%(2*q)<q){
      MPI_Send(Abuf+l+1, n-l-1, MPI_INT, myid+q, (i-1)*size+myid, MPI_COMM_WORLD);
      MPI_Probe(myid+q, (i-1)*size+myid+q, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_INT, &ncount);
      number_buf = (int*)malloc(sizeof(int) * (l+1+ncount));
      b = (int*)malloc(sizeof(int) * (ncount));
      MPI_Recv(b, ncount, MPI_INT, myid+q, (i-1)*size+myid+q, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      merge(l+1, ncount, Abuf, b, number_buf);
      free(Abuf);
      free(b);
      Abuf = number_buf;
      n = l+1+ncount;
    }else{
      MPI_Send(Abuf, l+1, MPI_INT, myid-q, (i-1)*size+myid, MPI_COMM_WORLD);
      MPI_Probe(myid-q, (i-1)*size+myid-q, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_INT, &ncount);
      number_buf = (int*)malloc(sizeof(int) * (ncount+n-(l+1)));
      b = (int*)malloc(sizeof(int) * (ncount));
      MPI_Recv(b, ncount, MPI_INT, myid-q, (i-1)*size+myid-q, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      merge(n-(l+1), ncount, Abuf+l+1, b, number_buf);
      free(Abuf);
      free(b);
      Abuf = number_buf;
      n = n-(l+1)+ncount;
    }   
  }

  return(0);
}
