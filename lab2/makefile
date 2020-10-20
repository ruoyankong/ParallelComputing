.SUFFIXES: .c .o .f .F

CC			=  mpicc  
CFLAGS			=  -g -Wall

FILES =  main.o aux.o hqsort.o

main.ex: $(FILES) 
	${CC} ${CFLAGS} -o main.ex $(FILES) -lm 

.c.o:
	${CC} ${CFLAGS} $< -c -o $@

clean:
	rm *.o *.ex
