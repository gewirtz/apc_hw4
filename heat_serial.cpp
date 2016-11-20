//Runs with ./heat_serial <nx> for a solution with grid size nx^2
//run until time t = 0.5 pi^2 / kappa
//2D domain of size 0 leq x leq pi and 0 leq y leq pi
//boundary conds: T(x,0) = cos^2x, T(x,pi)=sin^2x, T(0,y)=T(pi,y) (periodic in x)
//start from init cond T=0 everywhere. Run until t = 0.5pi^2/kappa
//deltat < deltax^2/4kappa

#include <iostream>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
using namespace std;

#define PI 3.14159265358979323846

int main(int argc, char* argv[]){

	struct timeval t_start,t_end;
	gettimeofday(&t_start, NULL);

	if(argc != 2){ //Check for the correct number of arguments
		cerr << "Usage: ./heat_serial <nx>" << endl;
		return(1);
	}
	const int nx = atoi(argv[1]);
	if(nx<1){ //check that the parameters given make sense
		cerr <<"Please enter positive nx " << endl;
		return(1);
	}


	//create the grids
	double **grid = new double*[nx]; 
	for(int i = 0; i < nx; i++){
		grid[i] = new double[nx];
	}

	double **nextgrid = new double*[nx];
	for(int i=0; i<nx; i++){
		nextgrid[i] = new double[nx];
	}
	
	//set initial and boundary conditions
	for(int i=0; i<nx; i++){
		for(int j=0; j<nx; j++){
			if(j==0){
				grid[i][j] = nextgrid[i][j] = cos(i*PI/double(nx)) * cos(i*PI/double(nx));
			}
			else if(j==(nx-1)){
				grid[i][j] = nextgrid[i][j] = sin(i*PI/double(nx)) * sin(i*PI/double(nx));
			}
			else{
				grid[i][j] = nextgrid[i][j] = 0;
			}
		}
	}

	//run the centered finite difference method
	double dx = PI/nx;
	const double kappa = 1;
	const double dt = (dx*dx)/(6*kappa);
	const double tlim = 0.5*(PI*PI)/kappa;
	const double stepnum = tlim/dt;
	double numerator;


	for(int k = 0; k < stepnum; k++){
		for(int i = 1; i < (nx-1); i++){
			for(int j = 1; j<(nx-1); j++){
				numerator = grid[i-1][j] + grid[i+1][j] + grid[i][j-1] + grid[i][j+1] -4*grid[i][j];
				nextgrid[i][j] = grid[i][j] + (dt*kappa*numerator/(dx*dx));
			}
		}
		
		for(int i=1; i<(nx-1);i++){ //for the first and last rows
			numerator = grid[0][i-1] + grid[0][i+1] + grid[nx-1][i] + grid[1][i] - 4*grid[0][i];
			nextgrid[0][i] = nextgrid[nx-1][i] = grid[0][i] + (dt*kappa*numerator/(dx*dx));
		}
		
		//change pointers from one array to the other
		for(int i=0; i<nx; i++){
			for(int j=0; j<nx; j++){
				grid[i][j] = nextgrid[i][j];
			} //this may not work...
		}	
	}

	//end the clock
	gettimeofday(&t_end,NULL);
	float elapsed = ((t_end.tv_sec -t_start.tv_sec) * 1000000u + t_end.tv_usec - t_start.tv_usec) / 1.e6;
	cout << "Time: " << elapsed << endl;
	//now get the final, volume-averaged temperature
	double total=0;
	for(int i=0;i<nx;i++){
		for(int j=0;j<nx;j++){
			total = total + grid[i][j];
		}
	}
	cout << "The final volume-averaged temperature is: " << total/(nx*nx) << endl;

	//generate file to draw plot from in R or something
	char fname[25];
	sprintf(fname,"serial_%d.txt",nx);
	ofstream fout(fname);
	for(int i=0; i<nx; i++){
		for(int j=0; j<nx; j++){
			fout << i <<" "<<j<<" "<<grid[i][j] << endl;
		}
	}
	fout << endl;
	fout.close();

	//free memory
	for(int i=0; i<nx; i++){
		delete [] grid[i];
		delete [] nextgrid[i];
	}
	delete [] grid;
	delete [] nextgrid;
	return(0);
}
