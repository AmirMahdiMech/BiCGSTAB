#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>
using namespace std;

#define File "e20r0100.mtx"

void sizeread(int &size, int &NZ) {

	std::ifstream fin(File);
	
	int M, N, L;
	
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	
	fin >> M >> N >> L;
	size = M;
	
	double* matrix;		    
	matrix = new double[M*N];	   
	std::fill(matrix, matrix + M*N, 0.); 
	

	for (int l = 0; l < L; l++)
	{
		int m, n;
		double data;
		fin >> m >> n >> data;
		matrix[(m-1) + (n-1)*M] = data;
	}
	int counter = 0;
	for (int m = 0; m < M; m++)
	{
	    for(int n = 0; n < N; n++){
	        if (matrix[m + n*M]!=0)
	        counter++;
	        if (m == n && matrix[m + n*M]==0){
	        matrix[m + n*M] = 1E-6;
	        counter++;}
	}
	}
	NZ = counter;
	fin.close();	
}

void Matrixread(double a[], int ra[], int ca[]){

	std::ifstream fin(File);
	
	int M, N, L;
	
	while (fin.peek() == '%') fin.ignore(2048, '\n');
	
	fin >> M >> N >> L;
	
	double* matrix;		     
	matrix = new double[M*N];	     
	std::fill(matrix, matrix + M*N, 0.); 
	
	for (int l = 0; l < L; l++)
	{
		int m, n;
		double data;
		fin >> m >> n >> data;
		matrix[(m-1) + (n-1)*M] = data;
	}
	int counter = 0;
	for (int m = 0; m < M; m++)
	{
	    for(int n = 0; n < N; n++){
	        if (matrix[m + n*M]!=0)
	        counter++;
	        if (m == n && matrix[m + n*M]==0){
	        matrix[m + n*M] = 1E-6;
	        counter++;}
	}
	}
	int c = 0;
    double threshold = 1E-7;
    ra[0] = 0;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            double val = matrix[m * N + n];
            if (val >= threshold && c < counter) {
                a[c] = val;
                ca[c] = n;
                c++;
            }
        }
        ra[m+1] = c;
    }	
	fin.close();	
}

void SparseAx(double a[],int ra[],int ca[],double x[], int N, double w[]){
	
	int i;
	int k;
	int k1;
	int k2;
 	for (int i=0;i<N;i++){
	    w[i]=0;             
	    k1 = ra[i];
	    k2 = ra[i+1];
	    for ( k = k1 ; k < k2 ; k++)
	    {
	    	w[i] = w[i] + a[k]*x[ca[k]];
		}
 	}
 		return;
}

double dotproduct(double a[],double b[],int N){
 	  double val=0;
      for (int i=0;i<N;i++){
          val=val+a[i]*b[i];
      }
 	return val;
}

void bicgstab(double a[],int ra[],int ca[],double rhs[],double x[], int N){

	double *r0;                  
	r0=new double [N];
	double *r;                   
	r=new double [N];
	double *r_old;                  
	r_old=new double [N];
	double *p;                  
	p=new double [N];
	double *s;                   
	s=new double [N];
	double *dummy;
	dummy=new double[N];
                
	for (int i=0;i<N;i++)
	{     
	    r0[i]=r[i]=r_old[i]=p[i]=s[i]=x[i]=0;
	}
	double alpha=0;             
	double omega=0;                  
	double beta=0;               

	SparseAx(a,ra,ca,x,N,dummy);
	for (int i=0;i<N;i++)
	{
	    r0[i]=rhs[i]-dummy[i];
	}

	for (int i=0;i<N;i++)
	{
	    p[i]=r0[i];
	    r[i]=r0[i];
	}

	int MAX_itr=0;
	double tol=0;
	cout<<"Enter the maximum number of Iterations =";
	cin>>MAX_itr;
	cout<<"\n"<<"Enter Tolerance ="; 
	cin>>tol;

	for (int itr=0;itr<MAX_itr;itr++){

    SparseAx(a,ra,ca,p,N,dummy);
    alpha=(dotproduct(r,r0,N))/(dotproduct(dummy,r0,N));
    
    for (int i=0;i<N;i++){
    s[i]=r[i]-alpha*dummy[i];
    }

    SparseAx(a,ra,ca,s,N,dummy);
    omega=(dotproduct(dummy,s,N))/(dotproduct(dummy,dummy,N));

    for (int i=0;i<N;i++){
    x[i]=x[i]+alpha*p[i]+omega*s[i];
    }

    for (int i=0;i<N;i++){
    r_old[i]=s[i]-omega*dummy[i];
    }

    beta=((dotproduct(r_old,r0,N))/(dotproduct(r,r0,N)))*(alpha/omega);

    SparseAx(a,ra,ca,p,N,dummy);
    for (int i=0;i<N;i++){
    p[i]=r_old[i]+beta*(p[i]-omega*dummy[i]);
    }

    for (int i=0;i<N;i++){
    r[i]=r_old[i];
    }

    double res=sqrt((dotproduct(r,r,N))/N);
    cout<<" Iteration: "<<itr+1<<"  Error: "<<res<<'\n';
    if (res<=tol){
    break;
    }
}

return;
}

void matrixToCSR(int &n, int &counter, double** matrix, int *ra, int *ca, double *a){
	
	counter = 0;
	cout<<"Your matrix is "<<n<<" by "<<n<<endl;
	cout<<"Now you have to input the elements: "<<endl;
	
    matrix = new double*[n];
    bool ok = false;

    while (!ok){
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout<<"a_"<<i+1<<j+1<<": ";
            cin >> matrix[i][j];
            cout<<endl;
            if(matrix[i][j]!=0){
            	counter++;
			}
			else if(i == j){
				matrix[i][j] = 1E-6;
			}
        }
    }
    cout<<"Your matrix is :"<<endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout<<matrix[i][j]<<"\t";
            
        }
            cout<<endl;
    }
    cout<<"If this is your matrix press 1, if not press 0 :"<<"\t";
    cin>>ok;
	}
    ra[0] = 0;
    int k = 0;
    for (int i = 0; i < n; i++) {

        ra[i] = k;
        for (int j = 0; j < n; j++) {
            if (matrix[i][j] != 0.0) {
                ca[k] = j;
                a[k] = matrix[i][j];
                k++;
            }
        }
    }
    ra[n] = k;  
    delete matrix;
}

void Yourtest(){
	int counter;
    int n;
	cout<<" Please enter the size of your matrix: ";
	cin>>n;
    double **matrix;
    int Max = n*n;
    int ra[n+1];
    double a[Max];
    int ca[Max];
    double b[n];
    double x[n];
    matrixToCSR(n,counter,matrix,ra,ca,a);
	cout<<" Input the RHS of the equation: "<<endl;
    for ( int i = 0; i < n; i++ )
    {
      cout<<"\t b["<<i+1<<"] = ";
      cin>>b[i];
      x[i] = 0.0;
    }
    bicgstab(a,ra,ca,b,x,n);
    cout<<" Final Solution: "<<endl;
    for ( int i = 0; i < n; i++ )
    {
      cout<<"\t x["<<i<<"] = "<<x[i]<<endl;
    }
}


void MTXtest(){
	
  int N , NZ_NUM;
  sizeread(N,NZ_NUM);
  double a[NZ_NUM];
  int ra[N+1];
  int ca[NZ_NUM];
  Matrixread(a, ra, ca);
  int itr_max;
  int i;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  double rhs[N];
  double x[N];
  double *x_estimate;
  x_estimate=new double [N];
  srand(time(NULL));
  for ( i = 0; i < n; i++ )
    {
      x[i] = 1.0;
      x_estimate [i]= rand()*1.0/RAND_MAX;
      
    }
  	SparseAx ( a, ra, ca, x, N, rhs);
  	int test;
  	double tol_abs;
  	double tol_rel;
  	double x_error;
  
  	cout << "\n";
  	cout << "  Test BICGSTAB on a sample matrix from matrixmarket.\n";
  	cout << "\n";

	bicgstab(a,ra,ca,rhs,x_estimate,N);


	for (i=0;i<N;i++){
	    cout<<"x["<<i+1<<"]= "<<x_estimate[i]<<"\n";
	    }

  	return ;
}


int main(int argc, char** argv)
{ 
    //Yourtest();
    MTXtest();
    return 0;
}





