//Runs with ./heat_serial <nx> for a solution with grid size nx^2
//run until time t = 0.5 pi^2 / kappa
//2D domain of size 0 leq x leq pi and 0 leq y leq pi
//boundary conds: T(x,0) = cos^2x, T(x,pi)=sin^2x, T(0,y)=T(pi,y) (periodic in x)
//start from init cond T=0 everywhere. Run until t = 0.5pi^2/kappa
//deltat < deltax^2/4kappa

#include <iostream>
#include <math>
using namespace std;

#define PI 3.14159265358979323846;

int main(int argc, char* argv[]){

	if(argc != 2){ //Check for the correct number of arguments
		cerr << "Usage: ./heat_serial <nx>" << endl;
		return(1);
	}
	const int nx = atoi(argv[1]);

	//create the grids
	double **grid = new double[nx]; 
	for(int i = 0; i < nx; i++){
		grid[i] = new double[nx];
	}

	double **nextgrid = new double[nx];
	for(int i=0; i<nx; i++){
		nextgrid[i] = new double[nx];
	}
	
	//set initial and boundary conditions
	for(int i=0; i<nx; i++){
		for(int j=0; j<nx; j++){
			if(j==0){
				grid[i][j] = nextgrid[i][j] = cos(i*PI/double(nx)) * cos(i*PI/double(nx));
			}
			else if(j==nx-1){
				grid[i][j] = nextgrid[i][j] = sin(i*PI/double(nx)) * sin(i*PI/double(nx));
			}
			else{
				grid[i][j] = nextgrid[i][j] = 0;
			}
		}
	}

	for(int i=0; i<nx; i++){
		grid[0][i] = nextgrid[0][i] = grid[nx-1][i];
	}

	//run the centered finite difference method
	double dx = PI/nx;
	const double kappa = 2;
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
		for(int i=1;i<(nx-1);i++){ //for the first column
			numerator = grid[i][nx-1] + grid[i][1] + grid[i-1][0] + grid[i+1][0] - 4*grid[i][j];
			nextgrid[i][0] = grid[i][j] + (dt*kappa*numerator/(dx*dx));
		}
		for(int i=1;i<(nx-1);i++){ // for the last column
			numerator = grid[i][nx-2] + grid[i][0] + grid[i-1][nx-1] + grid[i+1][nx-1] -4*grid[i][j];
			nextgrid[i][0] = grid[i][j] + (dt*kappa*numerator/(dx*dx));
		}
		for(int i=1; i<(nx-1);i++){ //for the first row
			numerator = grid[0][i-1] + grid[0][i+1] + grid[nx-1][i] + grid[1][i] - 4*grid[i][j];
			nextgrid[0][i] = grid[i][j] + (dt*kappa*numerator/(dx*dx));
		}
		numerator = grid[0][nx-1] + grid[0][1] + grid[nx-1][0] + grid[1][0] -4*grid[0][0]; //upper left corner
		nextgrid[0][0] = grid[0][0] + (dt*kappa*numerator/(dx*dx));
		numerator = grid[nx-2][0] + grid[0][0] + grid[nx-1][nx-1] + grid[nx-1][1] -4*grid[nx-1][0]; //lower left corner
		nextgrid[nx-1][0] = grid[nx-1][0] + (dt*kappa*numerator/(dx*dx));
		numerator = grid[0][nx-2] + grid[0][0] + grid[nx-1][nx-1] + grid[1][nx-1] -4*grid[0][nx-1]; //upper right corner
		nextgrid[0][nx-1] = grid[0][nx-1] + (dt*kappa*numerator/(dx*dx));
		numerator = grid[nx-1][nx-2] + grid[nx-1][0] + grid[nx-2][nx-1] + grid[0][nx-1] -4*grid[nx-1][nx-1]; //bottom right
		nextgrid[nx-1][nx-1] = grid[nx-1][nx-1] + (dt*kappa*numerator/(dx*dx));
		
		//change pointers from one array to the other
		for(int i=0; i<nx; i++){
			for(int j=0; j<nx; j++){
				grid[i][j] = nextgrid[i][j];
			} //this may not work...
		}	
	}

	//free memory
	for(int i=0; i<nx; i++){
		delete [] grid[i];
		delete [] nextgrid[i];
	}
	delete [] grid;
	delete [] nextgrid;
}
