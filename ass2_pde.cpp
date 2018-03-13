// Assignment 2
// Crank Nicolson Method
 
#include<iostream>
#include<math.h>
#include <cmath> 
using namespace std;

int N;
//#define N [10];

int main(void){
	
	N = 200; // number of sample point 
	double lamda = N;
	double h = 1.0/N; // grid size
	double tau = 1.0/N; // time grid size
	double T = 0.2;

	double TA[N-1][N-1],A[N-1][N-1],B[N-1][N-1]; //for TA for thomas algorithm, A = A_lambda, B = B_lambda
	int N_t = T/tau;
//	cout<<N_t;
	double M[N_t+1][N+1],M_exact[N_t+1][N+1]; 	// M is matrix for the whole field u, M_exact is the exact solution given by sin(pi*x)*exp(-pi^2*t)

// Boundry Condition
	for(int i=0;i<N_t+1;i++)
		{
		M[i][0] = 0.0;
		M[i][N] = 0.0;
		}
	for(int i=1;i<N;i++)
		{M[0][i] = sin(M_PI*i*h);}   // need to check if pi is included

// building the matrix A_lam total N-1 terms
	for (int i =0; i<N-1; i++){
		for (int j =0; j<N-1; j++){
			if (j == i-1){A[i][j] = -lamda/2.0;}
			else if (j == i){A[i][j] = 1.0 + lamda;}
			else if (j == i+1){A[i][j] = -lamda/2.0;}
			else{B[i][j] = 0.0;}
		}
	}

// building the matrix B_lam  // total N-1 terms
	for (int i =0; i<N-1; i++){
		for (int j =0; j<N-1; j++)
		{
			if (j == i-1){B[i][j] = lamda/2.0;}
			else if (j == i){B[i][j] = 1.0 - lamda;}
			else if (j == i+1){B[i][j] = lamda/2.0;}
			else{A[i][j] = 0.0;}
		}
	}

	double f[N-1]; // column matrix b for euation of type AX=b, here f= B(u^n)
	double X[N-1]; // solution for u^(n+1) for i = 1 to N-1
//--------------------------exact solution
	for (int i= 0;i<N_t+1;i++)
	{	for(int j = 0;j<N+1;j++)
		{	M_exact[i][j] = sin(M_PI*j*h)*exp(-pow(M_PI,2)*i*tau)	;

		}
	}
//---------------------------Numerical solution
for (int w = 1; w< N_t+1; w++)		// bigger outer for loop
{  								// matrix multiplication B_lamda and u^n
		for (int j = 0;j<N-1;j++)
			{
			double sum = 0;
			for (int k=0;k<N-1;k++)
				{
				sum += B[j][k]*M[w-1][k+1];
				}
			f[j] = sum;
			}

//----------------------------------------------------------
	// Thomas algorithm for solving tridiogonal matrix

		for (int i =0; i<N-1; i++)
		{
			for (int j =0; j<N-1; j++)
				{
				TA[i][j] = A[i][j];
				}	
		}						//now we have matrix TA = A

		int N2 = N-1;		// TA ---> N-1 element for convenience
// transformation in the matrix		
		for(int i =1;i<N2;i++){
			double temp = TA[i][i-1]/TA[i-1][i-1];
			TA[i][i] = TA[i][i] - temp*TA[i-1][i]; // formula_class notes
			f[i] = f[i] - temp*f[i-1]; //formula
                        TA[i][i-1]= 0.0;
		}
//-----------------------------------------------------------
	 //back-substitution
//		X[N2-1] = f[N2-1]/TA[N2-1][N2-1];
		for (int i=N2-1;i>=0;i--)       
		{	
			double temp = f[i] - TA[i][i+1]*X[i+1];
			X[i] = temp/TA[i][i];
		}

// filling X in M
		for (int l = 1;l<N;l++)
			{
			M[w][l] = X[l-1];
			}
} // end for loop

//-------------------------------------------------
		
	// M is the final matrix solution
	// printing the solution of in the matrix form

//error------------------- 
	double	temp = 0.0;
	double error;
	for (int i =0;i<N_t+1;i++)
		{ 
		for (int j =0;j<N+1;j++)
			{	error  = abs(M_exact[i][j]-M[i][j]);
				if (error > temp)
					{temp =  error;}
			}
		}

cout<< endl<<"the error is "<<temp<<endl<<endl;	
	
	
}

// end of program 
//sublime text
// color coding
