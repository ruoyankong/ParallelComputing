#include <stdio.h>
#include <sys/time.h>
#define MAX_LEN 524288
#define NITER 100
#define BSIZE 256

__global__ void saxpy_par(int vecLen, float a, float *x_d, float *y_d){ 
   int i = blockIdx.x * blockDim.x + threadIdx.x;]
   if (i < vecLen)
       y_d[i] = y_d[i] + a * x_d[i]; 
} 

float saxpy_check(int n, float a, float *x, float *y, float *z)
	// a, x, y == original data for saxpy
	// z = result found -- with which to compare.
	float s=0.0, t = 0.0;
	for (int i=0; i<n; i++) {
		y[i] += a * x[i] ;
		s += (y[i] - z[i])*(y[i] - z[i]);
		t += z[i]*z[i];
	}
	if (t == 0.0) return(-1);
	else
	return(sqrt(s/t));
}

int main (int argc, char *argv[]){
	float  *x, *y, *z;
	float a = 1.0;
	float error, t, mflops;
	int iter;
	x 	= (float*) malloc(MAX_LEN*sizeof(float));
	y 	= (float*) malloc(MAX_LEN*sizeof(float));
	z 	= (float*) malloc(MAX_LEN*sizeof(float));
	a = a/NITER;
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
		cudaMemcpy(xd, x, vecLen*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(yd, y, vecLen*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(ad, &a, sizeof(float), cudaMemcpyHostToDevice);
	    dim3 dimBlock(BSIZE);
  		dim3 dimGrid((vecLen + BSIZE-1) / BSIZE);
  		t = wctime();
  		for (iter = 0 ;iter<NITER; iter++)
			saxpy par <<<dimGrid, dimBlock>>> (vecLen, a, x_d, y_d);
		t = wctime() - t;
		cudaMemcpy(z, y_d, vecLen*sizeof(float), cudaMemcpyDeviceToHost);
		error = saxpy_check(vecLen, a, x, y, z);
		mflops = 2.0*vecLen/1000000.0/t;
		printf(" ** vecLen = %d, Mflops =%10.2e, err = %10.2e\n", vecLen, mflops, error); 
	}




}