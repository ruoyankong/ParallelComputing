###NVCC = /usr/local/cuda/bin/nvcc

NVCC = nvcc

##  NVCC_FLAGS = -I/usr/local/cuda/include
# make dbg=1 tells nvcc to add debugging symbols to the binary
##      NVCC_FLAGS += -g -O0

NVCC_FLAGS += -O3 -Wno-deprecated-gpu-targets
NVCC_FLAGS += -D_FORCE_INLINES

testSaxpy.ex:testSaxpy.cu
	$(NVCC) $(NVCC_FLAGS) -lcuda testSaxpy.cu -o testSaxpy.ex

testSaxpyC.ex:testSaxpyC.cu
	$(NVCC) $(NVCC_FLAGS) testSaxpyC.cu -o testSaxpyC.ex -lcuda

clean:
	rm *.ex
