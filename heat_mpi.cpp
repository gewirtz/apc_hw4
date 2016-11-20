//This program runs using MPI to parallelize the solution to the heat differential equation
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/time.h>

using namespace std;

#define PI 3.14159265358979323846


int main(int argc, char *argv[]){

	struct timeval t_start,t_end;
	gettimeofday(&t_start, NULL);

	if(argc != 2){ //Check for the correct number of arguments
		cerr << "Usage: ./heat_mpi <nx>" << endl;
		return(1);
	}
	const int nx = atoi(argv[1]);
	if(nx<1){ //check that the parameters given make sense
		cerr <<"Please enter positive nx " << endl;
		return(1);
	}

	int numtasks, rank, rc;
	rc=MPI_Init(&argc, &argv);
	if(rc!=MPI_SUCCESS){
		MPI_Abort(MPI_COMM_WORLD,rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//create the grids
	int nrows = (nx/numtasks)+ 2; //flip rows and columns
	double grid[nrows][nx]; 
	double nextgrid[nrows][nx];

	//set initial and boundary conditions
	for(int i=0; i<nrows; i++){
		for(int j=0; j<nx; j++){
			if((i==1) && (rank==0)){
				grid[i][j]  = cos(j*PI/double(nx)) * cos(j*PI/double(nx));
				nextgrid[i][j] = cos(j*PI/double(nx)) * cos(j*PI/double(nx));
			}
			else if((i==(nrows-2)) && (rank==(numtasks-1))){
				grid[i][j]  = sin(j*PI/double(nx)) * sin(j*PI/double(nx));
				nextgrid[i][j]= sin(j*PI/double(nx)) * sin(j*PI/double(nx));
			}
			else{
				grid[i][j]  = 0;
				nextgrid[i][j] = 0;
			}
		}
	}
	
	//run the centered finite difference method
	const double dx = PI/nx;
	const double kappa = 1;
	const double dt = (dx*dx)/(6*kappa);
	const double tlim = 0.5*(PI*PI)/kappa;
	const double stepnum = tlim/dt;

	double numerator;
	MPI_Request requests[8];
	MPI_Status stati[8];
	for(int k = 0; k < stepnum; k++){
		int start=1; //assume that we start updating on the second row (skip ghost row)
		int end=nrows-1; //to end before this so i < end
		if(rank ==0){
			start=2; //skip update of boundary conditions
		}
		if(rank==(numtasks-1)){
			end=nrows-2;
		}

		for(int i = start; i < end; i++){ //updates nextgrid all rows not ghost cells or boundary
			for(int j = 1; j<(nx-1); j++){ //skip update of first and last cols (no ghost cells)
				numerator = grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1] -4*grid[i][j];
				nextgrid[i][j] = grid[i][j] + (dt*kappa*numerator/(dx*dx));
			}
		}

		for(int i=start; i<end;i++){ //for the first and last columns (originally rows)
			numerator = grid[i-1][0] + grid[i+1][0] + grid[i][nx-1] + grid[i][1] - 4*grid[i][0];
			nextgrid[i][0] = nextgrid[i][nx-1] = grid[i][0] + (dt*kappa*numerator/(dx*dx));
		}

		for(int i=start; i<end; i++){ 
			for(int j=0; j<nx; j++){
				grid[i][j] = nextgrid[i][j]; //copy the nextgrid values to grid for non ghost/boundary cells
			}
		}

		int left = (rank-1+numtasks)%numtasks;
		int right = (rank+1)%numtasks;
		int p,p1,p2,p3;
		p=MPI_Isend(&grid[nx/numtasks][0],nx,MPI_DOUBLE,right,7,MPI_COMM_WORLD,&requests[0]);
		if(p!=MPI_SUCCESS){
			cerr << "Error in MPI_Isend" << endl;
			exit(1);
		}
		p1=MPI_Irecv(&grid[0][0],nx,MPI_DOUBLE,left,7,MPI_COMM_WORLD,&requests[1]);
		if(p1!=MPI_SUCCESS){
			cerr << "Error in MPI_Irec" << endl;
			exit(1);
		}
		p2=MPI_Isend(&grid[1][0],nx,MPI_DOUBLE,left,9,MPI_COMM_WORLD,&requests[2]);
		if(p2!=MPI_SUCCESS){
			cerr << "Error in MPI_Isend" << endl;
			exit(1);
		}
		p3=MPI_Irecv(&grid[nrows-1][0],nx,MPI_DOUBLE,right,9,MPI_COMM_WORLD,&requests[3]);
		if(p3!=MPI_SUCCESS){
			cerr << "Error in MPI_Irec" << endl;
			exit(1);
		}

		MPI_Waitall(4,requests,stati);

	}

	char fname[50];
	sprintf(fname,"mpi_%d_rank_%d.txt",nx,rank);
	ofstream fout(fname);
	double ptotal=0;
	for(int i=1; i<(nrows-1); i++){
		for(int j=0; j<nx; j++){
			fout << i <<" "<<j<<" "<<grid[i][j] << endl;
			ptotal = ptotal + grid[i][j];
		}
	}
	fout << endl;
	fout.close();

	ptotal = ptotal/((nrows-2)*nx);
	MPI_Request request[numtasks-1];
	//each process sends its average temp to main
	if(rank!=0){
		int p2;
		p2=MPI_Isend(&ptotal,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,&request[rank-1]);
		if(p2!=MPI_SUCCESS){
			cerr << "Error in MPI_Isend" << endl;
		return(1);
		}
	}
	MPI_Status status[numtasks-1];

	MPI_Barrier(MPI_COMM_WORLD);
	//rank==0 gathers them all and averages
	if(rank==0){
		int p3;
		double ptotals[numtasks];
		ptotals[0]=ptotal;
		for(int i=1;i<numtasks;i++){
			p3=MPI_Irecv(&ptotals[i],1,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&request[i-1]);
			if(p3!=MPI_SUCCESS){
				cerr << "Error in MPI_Irec" << endl;
				return(1);
			}
		}
		double ftotal=0;
		for(int i=0;i<numtasks;i++){
			ftotal=ftotal+ptotals[i];
		}
		cout << "The final volume-averaged temperature is: " << ftotal/(numtasks) << endl;
	}

	gettimeofday(&t_end,NULL);
	float elapsed = ((t_end.tv_sec -t_start.tv_sec) * 1000000u + t_end.tv_usec - t_start.tv_usec) / 1.e6;
	cout << "Time: " << elapsed <<" for rank: " << rank<< endl;
	MPI_Finalize();
	return 0;
}


