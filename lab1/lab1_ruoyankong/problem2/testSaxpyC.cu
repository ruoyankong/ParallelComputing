#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define MAX_LEN 524288
#define NITER 100
#define BSIZE 256


/*-------------------- POSIX-compliant timer in seconds */
double wctime() 
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec + tv.tv_usec/1000000.0);
}

__global__ void saxpy_par(int vecLen, float a, float *x_d, float *y_d){ 
   int i = blockIdx.x * blockDim.x + threadIdx.x;
   if (i < vecLen)
       y_d[i] = y_d[i] + a * x_d[i]; 
} 

float saxpy_check(int n, float a, float *x, float *y, float *z){
	// a, x, y == original data for saxpy
	// z = result found -- with which to compare.
	float s=0.0, t = 0.0;
	int i;
	for (i=0; i<n; i++) {
		y[i] += a * x[i] ;
		s += (y[i] - z[i])*(y[i] - z[i]);
		t += z[i]*z[i];
	}
	if (t == 0.0) return(-1);
	else
	return(sqrt(s/t));
}

int main (int argc, char *argv[]){
	float  *x, *y, *z, *xd, *yd, *ad;
	float a = 1.0;
	float error;
	long double mflops, t;
	int iter, i;
	x 	= (float*) malloc(MAX_LEN*sizeof(float));
	y 	= (float*) malloc(MAX_LEN*sizeof(float));
	z 	= (float*) malloc(MAX_LEN*sizeof(float));
	for (i=0; i<MAX_LEN; i++) {
	    x[i] = (float) rand() / (float) rand();
	    y[i] = (float) rand() / (float) rand();
	} 	

	if (cudaMalloc((void **) &xd, MAX_LEN*sizeof(float)) != cudaSuccess) 
	    exit(2);	      
	if (cudaMalloc((void **) &yd, MAX_LEN*sizeof(float)) != cudaSuccess) 
	    exit(3);	
	if (cudaMalloc((void **) &ad, sizeof(float)) != cudaSuccess) 
	    exit(4);

	for(int vecLen = 1024; vecLen <= MAX_LEN; vecLen*=2) {
		cudaMemcpy(ad, &a, sizeof(float), cudaMemcpyHostToDevice);
	    	dim3 dimBlock(BSIZE);
  		dim3 dimGrid((vecLen + BSIZE-1) / BSIZE);
  		t = wctime();
  		for (iter = 0 ;iter<NITER; iter++){
			cudaMemcpy(xd, x, vecLen*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(yd, y, vecLen*sizeof(float), cudaMemcpyHostToDevice);
			saxpy_par <<<dimGrid, dimBlock>>> (vecLen, a, xd, yd);
			cudaMemcpy(z, yd, vecLen*sizeof(float), cudaMemcpyDeviceToHost);
		}
		t = wctime() - t;
		error = saxpy_check(vecLen, a, x, y, z);
		mflops = 2.0*vecLen*iter/(1000000.0*t);
		printf(" ** vecLen = %d, Mflops =%10.2Le, err = %10.2e\n", vecLen, mflops, error); 
	}

	free(x); 
 	free(y);
 	free(z);
 	cudaFree(xd);
 	cudaFree(yd);


}
