#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

/* defalut size of matrix */
#define DEFAULT_SIZE 2000

/* LAPACK double precision matrix-matrix multiplication */
void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	    double *alpha,double *a, int *lda, double *b, int *ldb,
	    double *beta, double *c, int *ldc);

/* timing function */
double wctime();

/*--------------------call lapack to perform matrix_matrix_product 
    C(m by n) = A(m by k) * B(k by n)                   */

int mat_mat_product_lapack(int m, int n, int k, double *A, double *B,
			   double *C){
  /*-------------------- C = row-wise storage; Lapack=column-wise,
    Need to transpose in and out of Lapack */
  char yy = 'T';    
  double one = 1.0;
  double zero = 0.0;
  dgemm_(&yy, &yy, &m, &n, &k, &one, A, &m, B, &k, &zero, C, &m);
  return 0;
}

void transpose(double *B, double *B_COPY, int m, int n){
  int i, j;

  # pragma omp parallel for
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      B_COPY[j*n+i]=B[i*n+j];
    }
  }
  
}

/* mat_mat_product C(m by n) = A(m by k) * B(k by n)  */
int mat_mat_product(int m, int n, int k, double *A, double *B,
		    double *C){
  int nt, i1, i2, i3;
  #pragma omp parallel
  {
    nt = omp_get_num_threads();
  }
  double *B_T;
  B_T = (double*) malloc(k*n*sizeof(double));
  transpose(B, B_T, k, n);
  printf("Run with %d threads\n",nt);
  /* loop to perform:   C(i,j) = sum_t A(i,t) * B(t,j)
   * Outer loop i1: rows of A
   */
  /*-------------------- set all of C to zero first*/

  /*-------------------- Inner loop, linear comb of rows  of B */
  # pragma omp parallel
  for( i1 = 0 ; i1 < m*n ; i1 ++) {
    C[i1] = 0.0;
  }
  /*--------------------row i1 of C == lin comb. of rows of B*/

  # pragma omp parallel
  for( i1 = 0 ; i1 < m ; i1 ++)	  {
      for( i3 = 0 ; i3 < n ; i3 ++ ) {
	  for( i2 = 0 ; i2 < k ; i2 ++) {
	     C[i3+i1*n] += A[i2+i1*k]*B_T[i2 +i3*n];
	  }
      }
   }
  return 0;
}

/*-------------------- driver                     */
int main (int argc, char *argv[]){
  omp_set_num_threads(32);
  int	  i, j;
  int     m, n, k; 		//problem sizes
  double  *a, *b, *c, *c1;   	//matrices
  double  t1, value, err;       // for checking errors
  
  /*-------------------- set matrix dimension */
  m = n = k = DEFAULT_SIZE;
  
  /*-------------------- set initial values */
  a 	= (double*) malloc(m*k*sizeof(double));
  b 	= (double*) malloc(k*n*sizeof(double));
  c 	= (double*) malloc(m*n*sizeof(double));
  c1	= (double*) malloc(m*n*sizeof(double));
  
  /*-------------------- note: LAPACK is FORTRAN based - assumes
    matrix stored by columns -- Lapack must use transposition.
    set A  */
  for(i = 0 ; i < m ; i ++)	{
    for(j = 0 ; j < k ; j ++) {
      value = (double)rand() / (double)RAND_MAX;
      a[j+i*k]  = value;
    }
  }
  
  /*-------------------- set B */
  for(i = 0 ; i < k ; i ++) {
    for(j = 0 ; j < n ; j ++)  {
      value = (double)rand() / (double)RAND_MAX;
      b[j+i*n] 	= value;
    }
  }
  
  /*-------------------- call Lapack - get the LAPACK time */
  t1 = wctime();
  mat_mat_product_lapack(m, n, k, a, b, c1);
  t1 = wctime()-t1;
  printf(" Lapack dgemm elapsed time %e\n",t1);

  /*-------------------- get time for your mat-mat function */
  t1 = wctime();
  mat_mat_product(m, n, k, a, b, c);
  t1 = wctime()-t1;

  /*-------------------- check error. */
  err = 0.0;
  for(i = 0 ; i < m ; i ++){
    for(j = 0 ; j < n ; j ++)	{
      err += (c1[i+j*m]-c[j+i*n])*(c1[i+j*m]-c[j+i*n]);
    }
  }
  
  printf(" My Matrix-matrix elapsed time %e, err is %e\n",t1,err);
  double flop_rate = (double)(2*m/100.0*n/100.0*k/100.0+m*n/1000000.0+k*n/1000000.0)/(t1);
  printf(" My flop rate is %e MFLOPS\n", flop_rate);
  
  /*-------------------- free arrays */
  free(a);
  free(b);
  free(c);
  free(c1);	
  return 0;
}
