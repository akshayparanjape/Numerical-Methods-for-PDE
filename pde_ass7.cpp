// Assignment 7 Numerical Methods PDE 
// Ritz Galerkin Method

#include <iostream>
#include <cmath>
using namespace std;

// solving -u^2 = f

#define N 1000

double x[N+1],u[N+1],u_exact[N+1];
double A[N-1][N-1],f[N-1],xsol[N-1];

//function on the right side f
double f_right(double x1){
	return M_PI*M_PI*sin(M_PI*x1); // check
}

double trap(double x[], int i){
	double part1 = (x[i] - x[i-1])/2.0*(f_right(x[i]) + f_right(x[i-1]));
	double part2 = (x[i+1] - x[i])/2.0*(f_right(x[i+1]) + f_right(x[i]));
	return 0.5*(part1 + part2);
}

int main(){

	// bondry condition
	u[0] = 0.0;u[N] = 0.0;

	double alpha = 10.0;
	x[0] = 0.0;
	x[N]=1.0;

	for (int i = 1; i < N; ++i)
	{	
		double t = N;
	 	//x[i] = i/t; // equal interval grids
	 	x[i] = (exp(alpha*i/N)-1.0)/(exp(alpha)-1.0); 
	}
// forming the tridiagonal matrix
	for (int j = 0; j < N-1; ++j)
	{	int i = j+1;
	 	A[j][j] = 1.0/(x[i+1] - x[i]) + 1.0/(x[i] - x[i-1]) ;
	 	if (j != N-1 ){A[j][j+1] = -1.0/(x[i+1] - x[i]);}
	 	if (j != 0){A[j][j-1] = -1.0/(x[i] - x[i-1]); }
	} 
//integration for the function on the right side
	// integration is done using trapezoidal rule
	for (int i = 0; i < N-1; ++i)
	{
		f[i] = trap(x,i+1);
	}

//---------Thomas algorithm---------------used in assignment #2

	int n = N-1;		
// transformation in the matrix		
		for(int i = 1;i<n;i++){
			double temp = A[i][i-1]/A[i-1][i-1];
			A[i][i] = A[i][i] - temp*A[i-1][i]; // formula_class notes
			f[i] = f[i] - temp*f[i-1]; //formula
                        A[i][i-1]= 0.0;
		}
//-----------------------------------------------------------
	 //back-substitution
//		X[n-1] = f[n-1]/TA[n-1][n-1];
		for (int i=n-1;i>=0;i--)       
		{	
			double temp = f[i] - A[i][i+1]*xsol[i+1];
			xsol[i] = temp/A[i][i];
		}

// updating the solution
	for (int i=1;i<N;i++){
		u[i] = xsol[i-1];
	}

// ----------------exact solution & error
	double error[N+1];
	double error_norm2 = 0.0;

	for (int i = 0; i < N+1; ++i)
	{
		u_exact[i] = sin(M_PI*x[i]);
		error[i] = abs(u[i] - u_exact[i]);
		if (error[i] > error_norm2){error_norm2 = error[i];}
	}

//----------printing-------------------	

	cout<<"Maximum Error "<<endl;
	cout<<error_norm2<<endl;
	double h[N];
	double h_max =0.0;

	for (int i = 0; i < N+1; ++i)
	{
		h[i] = x[i+1] - x[i];
		if (h[i] > h_max) {h_max = h[i];}
	}
	cout<<"h_max : \t "
	cout<<h_max<<endl;
	cout<<endl;

}
