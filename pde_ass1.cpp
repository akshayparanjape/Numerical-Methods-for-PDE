# include<iostream>
#include <cmath> 
using namespace std;

int p,q,m;

//power function
double pow(int x, int n){
	if (n==0){
		return 1;
	}else{
		return x*pow(x,n-1);
	}		
}

// factorial function
double fact(int n){
	if (n == 0){
		return 1;
	}else{
		return n*fact(n-1);
	}
}


int main(){
//	p=0;q=2;m=1;	

	cout<<"Enter the value for p"<<endl;
	cin>>p;
	cout<<"Enter the value for q"<<endl;
	cin>>q;
	cout<<"Enter the value for m"<<endl;
	cin>>m;

// check for the sufficiently many points
	if (p+q < m){
		cout<<"Do not have Sufficiently many point to support the approximation"<<endl;}
	else{

// the main program
	int T = p+q+1;  // T gives total number of rows of taylor table
	double a[T][T+1], x[T];	

// creating taylor table

// a is taylor tabble

	for(int i=0; i <T; i++){
		for (int j=0;j<T;j++){
			a[j][i] = pow(i-p,j)/fact(j);  // the element here are the coefficient for variable a,b,c
		}
	}
	// here we will get the transpose of actual taylor table, which will make it easy for gauss elimination

// defining matrix b in | Ax=b
	for (int i=0;i<T;i++){	a[i][T] = 0;}  

	a[m][T] = 1;
	

	int n=T;

// now solving using gthe gauss elimination method
	for (int i=0;i<n;i++)                    //Pivotisation
        for (int k=i+1;k<n;k++)
            if (a[i][i]<a[k][i])
                for (int j=0;j<=n;j++)
                {
                    double temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
//----------------------------------------	
	for (int i=0;i<n-1;i++)            //loop to perform the gauss elimination
        for (int k=i+1;k<n;k++)
            {
                double t=a[k][i]/a[i][i];
                for (int j=0;j<=n;j++)
                    a[k][j]=a[k][j]-t*a[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
            }
//-----------------------------------------------    

    for (int i=n-1;i>=0;i--)                //back-substitution
    {                        //x is an array whose values correspond to the values of a,b,c,d,e..
        x[i]=a[i][n];              
        for (int j=0;j<n;j++)
            if (j!=i)            
                x[i]=x[i]-a[i][j]*x[j];
        x[i]=x[i]/a[i][i];            //now finally divide the RHS by the coefficient of the variable to be calculated
    }
//---------------------------------------------------
    cout<<"\nThe values of the variables are as follows:\n";
    for (int i=0;i<n;i++)
        cout<<x[i]<<endl;            // Print the values of a,b,c,d,e....    
	cout<<endl;
// end of gauss elimination


// finding the order

	double order_mat[T][2]; // extra two rows in taylor table, since the order will be among one of this row

	for(int i =0;i<T;i++){
		order_mat[i][0] = pow(i-p,T)/fact(T);
		order_mat[i][1] = pow(i-p,T+1)/fact(T+1);
	}
// sum 1, sum2 are defined to check if the sumation of the coefficent with variables goes to zero

	double sum1 = 0;	double sum2=0;

	for (int i=0;i<T;i++){
		sum1 = sum1 + order_mat[i][0]*x[i]; // T+1 coloum in taylor table

		sum2 = sum2 + order_mat[i][1]*x[i]; // T+2 coloum in taylor table
	}

	int r; // order

	if (abs(sum1) > 1e-14){			
			r=T-m;
			cout<<"accuracy of the approximation, value of r_order  " <<r<<endl;	
		}
	else if (abs(sum2) > 1e-14){
			r=T-m+1;
			cout<<"accuracy of the approximation, value of r_order " <<r<<endl;	
		}
	else{
		r = 0;  // this imples all the subsequesnt term gives zero summation No order of error or no error
		cout<<"accuracy of the approximation, value of r_order " <<r<<endl;}
	}
}
