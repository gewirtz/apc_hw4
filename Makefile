#OBJS: heat_serial heat_omp heat_mpi
#targets: heat_serial.cpp heat_omp.cpp

CXX = g++
CXXFLAGS = -Wall  


all : serial omp mpi

serial : heat_serial.cpp 
	$(CXX) $(CXXFLAGS) heat_serial.cpp -o heat_serial

omp : heat_omp.cpp
	$(CXX) -fopenmp $(CXXFLAGS) heat_omp.cpp -o heat_omp

mpi : heat_mpi.cpp
	mpic++ heat_mpi.cpp -o heat_mpi

clean:
	rm heat_serial heat_omp heat_mpi





