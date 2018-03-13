// Assignment 4 PDE
// Upwind Scheme

#include<iostream>
using namespace std;
#include<cmath>
#include <fstream> // to get the output in a file

#define N 10

float h, tau, v,a,s,s2;

int M; //n*tau = T
float T;

int main(){

	v =0.9;
	h = 1.0/N;
	a = 1.0;
	tau = v*h/a;
	M = N/v;

// Matrix for storing 																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																						
	float M_upward[M+1][N];

//	float M_lax[M+1][N];

//Upwind scheme
	s = abs(v);
//Lax-Wendroff scheme
	s2 = pow(v,2);

//initial condition
	for (int i=0;i<N;i++)
	{
		M_upward[0][i] = sin(2*M_PI*i*h);

		//M_lax[0][i] = sin(2*M_PI*i*h);
	}

//Discontinous

/*
	for (int i=0;i<N/2;i++)
	{
	//	M_upward[0][i] = 1.0;
		M_lax[0][i] = 1.0;
	}
	for(int i=N/2;i<N;i++)
	{
	//	M_upward[0][i] = 0.0;
		M_lax[0][i] = 0.0;	
	}
*/

	for (int j = 1;j<M+1;j++)
	{

		for (int i = 0; i < N; i++)
		{
			int a = i-1;int b = i+1;

			if (a < 0){a = N + a;}
			else if (a>N-1){ a = a - N;}

			if (b < 0){b = N + b;}
			else if (b>N-1){ b = b - N;}

			M_upward[j][i] = (s+v)/2.0*M_upward[j-1][a] + (1.0-s)*M_upward[j-1][i] + (s-v)/2.0*M_upward[j-1][b];
	//		M_lax[j][i] = (s2+v)/2.0*M_lax[j-1][a] + (1.0-s2)*M_lax[j-1][i] + (s2-v)/2.0*M_lax[j-1][b];
		
		}
	}



//modifield upwind Scheme
float meu = abs(a)*h/2.0*(1-abs(v));

// error 
	float error = 0.0,error2 = 0.0,temp,temp2;
	
		
			for (int i = 0; i < N; i++)
			{
				temp = abs(M_upward[M][i]-sin(2.0*M_PI*(i*h - M*tau)));
//				temp2 = abs(M_lax[M][i]-sin(2.0*M_PI*(i*h - M*tau)));
//				float exact = sin(2*M_PI*(i*h - M*tau))*exp(-meu*4*pow(M_PI,2)*M*tau);
//				temp = abs(M_upward[M][i]-exact);
				if (temp > error){error = temp;}
		//		if (temp2 > error2){error2 = temp2;}
			}
		
	cout<<"The Upwind Scheme error for  "<< N<< "  is  "<<error<<endl;
	//cout<<"The Lax-Wendroff Scheme error for  "<< N<< "  is  "<<error2<<endl;

// Savinf the data in outpur file
	/*
ofstream out ("output_1000.csv");
for (int i = 0; i < N; ++i)
{
	out<< M_lax[0][i]<<","<<M_lax[M][i]<<endl;
}

out.close();
	*/
}
