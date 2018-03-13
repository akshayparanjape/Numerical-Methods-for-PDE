// Assignment No. 5
// Burger Equation

#include<iostream>
using namespace std;
#include<cmath>
#include <fstream> // to get the output in a file

#define N 100
#define epsilon 0.5e-5
#define T 1.0 // maximum permisiible time


float h, tau,t;//t denotes time, tau is time interval
float u[N],u_sol[N]; // u initial u_sol final solution 

float f( float u[N],int i){
	return 1.0/2.0*u[i]*u[i];


}

int main(){		

	h = 1.0/N;
	float t = 0.0; //time t = t + tau
	

// stop the iteration when tau >= T	

//initial condition (1)

	for (int i = 0; i < N; ++i)
		{
			u[i] = 1.0 + 1.0/2.0*sin(2.0*M_PI*(i*h));
		}	
		
//initial condition (2)
/*
	for (int i = 0; i < N; ++i)
		{
			if (i*h>0.1 && i*h<0.3){u[i] = 1;}
			else {u[i] = 0.0;}
		}
*/

//------------------------------------

	float a_arr[N]; // a_arr- = a_arr(i-1/2) || a_arr+ = a_arr(i+1/2)
	
	while (t < T) 
	{	
		for (int i = 0; i < N; ++i)
		{
			int a = i+1;
			if (a>N-1){ a = a - N;}

			if (abs(u[i]-u[a]) < epsilon)
				a_arr[i] = u[i]; //a_arr[i] = a_arr[i + 1/2]
			else
				a_arr[i] = (f(u,a)- f(u,i))/(u[a]-u[i]);
		
		}
		// CASE 1 alpha  = abs(a_arr)
		// calculating c+, c- for getting time interval tau
		float cplus,cminus;
		float c_max = 0.0;
		for (int i = 0; i < N; ++i)
		{
			cplus = abs(a_arr[i]) - a_arr[i];//cplus = c(i,i+1)  cminus = c(i,i-1)

			int m = i-1;
			if (m < 0){m = N + m;}
			cminus = abs(a_arr[m]) + a_arr[m];	
			float c = cplus + cminus;
			if (c > c_max){c_max = c;}
		} //c matrix done


		tau = h/c_max; // time interval
		// here I am taking the limiting condition for the value of tau
		
		float f_vector[N];
		for (int i = 0; i < N; ++i)
		{	int a = i+1;
			if (a>N-1){ a = a - N;}

			f_vector[i] = 0.5*(f(u,a) + f(u,i)) - abs(a_arr[i])/2.0*(u[a] - u[i]); //f_vecot0r[i] = i+1/2
		}

		for (int i = 0; i < N; ++i)
		{	int m = i-1;
			if (m < 0){m = N + m;}
			u_sol[i] = u[i] - tau/h*(f_vector[i] - f_vector[m]);
		}

	// modifying the array u so as to use in the next iteration loop
		for (int i = 0; i < N; ++i){u[i] = u_sol[i];}


		t = t + tau;

	}

// exporting the output in .csv file 
	// plotting has been done in python
ofstream out ("output.csv");
	for (int i = 0; i < N; ++i)
	{
		out<<u_sol[i]<<"  "<<endl;
	}
out.close();

}


