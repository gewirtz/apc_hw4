objects: heat_serial heat_omp
targets: heat_serial.cpp heat_omp.cpp

CC = g++
CFLAGS = -Wall 
O = -fopenmp

serial : heat_serial.cpp 
	$(CC) $(CFLAGS) $(^) -o $(@)

omp : heat_omp.cpp
	$(CC) -fopenmp $(CFLAGS) $O heat_omp.cpp -o heat_omp

clean:
	\rm heat_serial heat_omp



