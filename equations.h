#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
//#include <iostream>
//#include <malloc.h>
#include <time.h>
#include <math.h>
#include "/usr/include/complex.h"
#include <fftw3.h>

const double pi=3.141592653589793;

using namespace std;



    
void lintimes(double *times, double tmin, double tmax, double dt, double* delta_t, int nstat)
{
	times[0] = tmin;
	delta_t[0] = tmin/dt;
	double step = double(tmax)/double(nstat);
	for (int i=1;i<nstat;i++)
	{
		times[i] = tmin + i*step;
		delta_t[i] = (times[i]-times[i-1])/dt;
	}
	
	return;
}

//Writes a vector of logarithm-spaced times from tmin to tmax with nvolte points, and such that the difference beetween two consecutive times is a multiple of dt
void logtimes(double *times, double tmin, double tmax, double dt, double *delta_t, int nstat)
{
    double dn_stat;
    dn_stat = log10(tmax/tmin)/nstat;
    //printf("%le \n",dn_stat);
    if(tmin < dt)
    {
    	times[0] = tmin;
    	delta_t[0] = 1.0;
    }
    else
    {
    	times[0] = ceil(tmin/dt)*dt;
    	delta_t[0] = (times[0])/dt;
    }
    //delta_t[0] = 1.0;
    //times[0] = tmin;
    int i=1;
    
    for(i=1;i<nstat;i++)
    {
        //delta_t[i] = ceil( times[i-1]*(pow(10.0, dn_stat)-1)/dt );
        //times[i] = times[i-1] + dt*delta_t[i];
        times[i] = times[0]*pow(10, dn_stat*(i+1));
        if (times[i]-times[i-1] < dt)
        {
            times[i] = times[i-1] + dt;
            delta_t[i] = 1;
        }
        else
        {
            
            times[i] = ceil(times[i]/dt)*dt;
            delta_t[i] = (times[i]-times[i-1])/dt;
        }
        //printf("%le  %le \n",times[i], delta_t[i]);
       
    }
    return;
}

//Copy vector v1 into vector v2
void copy_vec(double* v1, double* v2, int sz)
{
    for(int i=0;i<sz;i++)
    {
        v2[i] = v1[i];
    }
    
    return;
}

void copy_vec(fftw_complex* v1, fftw_complex* v2, int sz)
{
    for(int i=0;i<sz;i++)
    {
        v2[i] = v1[i];
    }
    
    return;
}

//Copy one matrix into another
void copy_mat2mat(double* M1, double* M, int rgstart, int rgend, int colstart, int colend, int sz)
{
	int sz1 = colend-colstart+1;
	for(int rg=rgstart;rg <= rgend; rg++)
	{
		for(int col=colstart; col <= colend; col++)
		{
			M1[(rg-rgstart)*sz1 + (col-colstart)] = M[rg*sz+col];
		}
	}
	
	return;
}


//Copy vector v1 in the rg-th row of the matrix M
extern void copy_vec2matr(double* v1, double* M, int rg, int start, int end)
{
	int sz = end-start+1;
    for(int i=start;i<=end;i++)
    {
        M[rg*sz + i-start] = v1[i];
    }
    
    return;
}

//Copy vector v1 in the rg-th row of the matrix M
void copy_vec2matr(fftw_complex* v1, fftw_complex* M, int rg, int sz)
{
    for(int i=0;i<sz;i++)
    {
        M[rg*sz + i] = v1[i];
    }
    
    return;
}

//Extract row from a matrix and put it in a vector
void extract_row(double* v1,double* M, int rg, int colstart, int colend, int sz)
{
    //double *v1;
    //v1 = (double*) malloc(ncol*sizeof(double));
    for(int i=colstart;i<=colend;i++)
    {
        v1[i] = M[rg*sz + i];
    }
    
    return;
    
}

//Print vector to stdoutput
void print_vec(double* v1, int sz)
{
    for(int i=0;i<sz;i++)
    {
        printf("%le ",v1[i]);
    }
    printf("\n");
    return;
}

//Print vector to stdoutput
void print_vec(fftw_complex* v1, int sz)
{
    for(int i=0;i<sz;i++)
    {
        printf("%e + j*%e     \n",creal(v1[i]),cimag(v1[i]));
    }
    printf("\n");
    return;
}

void print_matr(double* M, int nrg, int ncol)
{
    for(int i=0;i<nrg;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            printf("%le ",M[i*ncol + j]);
        }
        printf("\n");
    }
    
    return;
}

void print_matr(fftw_complex* M, int nrg, int ncol)
{
    for(int i=0;i<nrg;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            printf("%le + i%le ",creal(M[i*ncol + j]),cimag(M[i*ncol + j]));
        }
        printf("\n");
    }
    
    return;
}
//Sums elements of a vector
double sum(double* v1, int start, int end)
{
    double S = 0.;
    int sz = end-start+1;
    for (int i=start;i<=end;i++)
    {
        S += v1[i];
    }
    
    return S;
}


//Mean of the elements of a vector
double mean(double* v1, int start, int end)
{
    double m = 0.;
    int sz = end-start+1;
    m = sum(v1, start, end)/sz;
    
    return m;
}

//Mean over the rows of a matrix
void mean_matrows(double* v, double* M, int nrg, int ncol)
{
	double m = 0.0;
	for (int col=0; col< ncol; col++)
	{
		for (int rg=0; rg< nrg; rg++)
		{
			m += M[rg*ncol + col];
		}
		v[col] = (double)m/nrg;
		m = 0.0;
	}
	
	return;
}

//Adds a scalar to every element of a vector
void scadd_vec(double* v1, double* v2, double x, int start, int end)
{
	for(int i=start;i<=end;i++)
	{
		v1[i] = v2[i] + x;
	}
	
	return;
}

//Multiply a vector by a scalar
void scmult_vec(double* v1, double* v2, double x, int start, int end)
{
	for(int i=start;i<=end;i++)
	{
		v1[i] = x*v2[i];
	}
	
	return;
}

void scmult_vec(fftw_complex* v1, fftw_complex* v2, double x, int start, int end)
{
	for(int i=start;i<=end;i++)
	{
		v1[i] = x*v2[i];
	}
	
	return;
}

void zeros(double *v1, int sz)
{
	for (int i=0;i<sz;i++)
	{
		v1[i] = 0.0;
	}
	
	return;
}

void zeros(fftw_complex *v1, int sz)
{
	for (int i=0;i<sz;i++)
	{
		v1[i] = 0.0 + I*0.0;
	}
	
	return;
}

void zeros(fftw_complex *v1, int start, int end)
{
	for (int i=start;i<=end;i++)
	{
		v1[i] = 0.0 + I*0.0;
	}
	
	return;
}

void zeros(double *v1, int start, int end)
{
	for (int i=start;i<=end;i++)
	{
		v1[i] = 0.0 ;
	}
	
	return;
}

void check(double *v, int i, double t)
{
	if(isnan(v[i]) || (!finite(v[i])))
	{
		printf("NaN encountered during integration: index %d, time %e \n", i, t);
		fflush(stdout);
		exit(1);
	}
	
	return;
}

void check(double v, int i, double t)
{
	if(isnan(v) || (!finite(v)))
	{
		printf("NaN encountered during integration: index %d, time %e \n", i, t);
		fflush(stdout);
		exit(1);
	}
	
	return;
}

void check(double v)
{
	if(isnan(v) || (!finite(v)))
	{
		printf("NaN encountered during integration\n ");
		fflush(stdout);
		exit(1);
	}
	
	return;
}



void print_vec2file(double* v1, double* v2, int sz, FILE* fout)
{
    double x, y;
    for(int i=0; i<sz;i++)
    {
        x = v1[i];
        y = v2[i];
        fprintf(fout, "%le %le \n", x, y);
    }
    
    return;
}

void print_vec2file(double* v1, int start, int end, FILE* fout)
{
    double x;
    for(int i=start; i<=end;i++)
    {
        x = v1[i];
        fprintf(fout, "%e ", x);
    }
    fprintf(fout, "\n");
    fflush(fout);
    
    return;
}

void print_vec2file_column(double* v1, int start, int end, FILE* fout)
{
    double x;
    for(int i=start; i<=end;i++)
    {
        x = v1[i];
        fprintf(fout, "%e \n", x);
    }
    fprintf(fout, "\n");
    fflush(fout);
    
    return;
}

void print_vec2file(fftw_complex* v1, int start, int end, FILE* fout)
{
    double x;
    double y;
    for(int i=start; i<=end;i++)
    {
        x = creal(v1[i]);
        y = cimag(v1[i]);
        fprintf(fout, "%e+%ei  ", x, y);
    }
    fprintf(fout, "\n");
    fflush(fout);
    
    return;
}

void print_mat2file(double* M, int nrg, int ncol, FILE* fout)
{
	double x;
	check(M[1]);
	for(int rg=0;rg<nrg;rg++)
	{
		for(int col=0;col<ncol;col++)
		{
			x = M[rg*ncol + col];
			fprintf(fout,"%le ",x);
			
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n");
	fprintf(fout,"\n");
	fflush(fout);
	
	return;
}

void print_mat2file_transpose(double* M, int nrg, int ncol, FILE* fout)
{
	double x;
	check(M[1]);
	for(int rg=0;rg<nrg;rg++)
	{
		for(int col=0;col<ncol;col++)
		{
			x = M[col*ncol + rg];
			
			fprintf(fout,"%le ",x);
			
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n");
	fprintf(fout,"\n");
	fflush(fout);
	
	return;
}




void print_mat2file(double* M, int rgstart, int rgend, int colstart, int colend, int sz, FILE* fout)
{
	double x;
	check(M[1]);
	for(int rg=rgstart;rg <= rgend;rg++)
	{
		for(int col=colstart;col <= colend;col++)
		{
			x = M[rg*sz + col];
			fprintf(fout,"%le ",x);
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n");
	fprintf(fout,"\n");
	fflush(fout);
	
	return;
}


void print_mat2file(fftw_complex* M, int nrg, int ncol, FILE* fout)
{
	double x,y;
	check(creal(M[1]));
	for(int rg=0;rg<nrg;rg++)
	{
		for(int col=0;col<ncol;col++)
		{
			x = creal(M[rg*ncol + col]);
			y = cimag(M[rg*ncol + col]);
			fprintf(fout,"%le %le  ",x,y);
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n");
	
	return;
}










double laplx(double* v, int rg, int col, int sz, double nux)
{
	double x = 0.0;
	x = nux*(v[rg*sz + col+1] + v[rg*sz + col-1] - 2.*v[rg*sz + col]);
	return x;
}

double laply(double* v, int rg, int col, int sz, double nuy)
{
	double x = 0.0;
	x = nuy*(v[(rg+1)*sz + col] +  v[(rg-1)*sz + col] - 2.*v[rg*sz + col]);
	return x;
}

double lapl2x(double* v, int rg, int col, int sz, double kx)
{
	double x = 0.0;
	x = -kx*(6*v[rg*sz + col] - 4*( v[rg*sz + col+1] + v[rg*sz + col-1] ) + v[rg*sz + col+2] + v[rg*sz + col-2]);
	return x;
}

double lapl2y(double* v, int rg, int col, int sz, double ky)
{
	double x = 0.0;
	x =  -ky*(6*v[rg*sz + col] - 4*( v[(rg+1)*sz + col] + v[(rg-1)*sz + col] ) + v[(rg+2)*sz + col] + v[(rg-2)*sz + col]);
	return x;
}

double lapl2xy(double* v, int rg, int col, int sz, double kxy)
{
	double x = 0.0;
	double a1 = 0.0;
	double a2 = 0.0;
	double a3 = 0.0;
	double a4 = 0.0;
	a1 = v[(rg+2)*sz + col+1] - v[(rg+2)*sz + col];
	a2 = v[(rg+1)*sz + col + 2] - v[(rg+1)*sz + col-1] + 6*(v[(rg+1)*sz + col] - v[(rg+1)*sz + col +2]);
	a3 = v[rg*sz + col -1] - v[rg*sz + col + 2] +  6*(v[rg*sz + col+1] - v[rg*sz + col]);
	a4 = v[(rg-1)*sz + col] - v[(rg-1)*sz + col+1];
	x = -kxy*(a1 + a2 + a3 + a4);
	
	return x;
}

double f(double* v, int rg, int col, int sz, double dx, double dt, double nux, double nuy, double kx, double ky, double kxy, double D, double noise)
{
	double x = 0.0;
	double d2x = 1.0/dx/dx;
	double d4x = 1.0/dx/dx/dx/dx;
	double C3 = sqrt(D*dt/dx);
	
	x = d2x*(laplx(v,rg,col,sz,nux) + laply(v,rg,col,sz,nuy)) + d4x*(lapl2x(v,rg,col,sz,kx) + lapl2y(v,rg,col,sz,ky) + lapl2xy(v,rg,col,sz,kxy)) + C3*noise;
	
	return x;
}

//Variance of the elements of a vector
double var(double* v1, int sz, double dx)
{
	//print_vec(v1, sz);printf("\n");
	int L = int(sz*dx);
    double m = mean(v1, 0, sz-1);
    double v = 0.0;
    for(int i=0;i< sz;i++)
    {
        v += (v1[i]-m)*(v1[i]-m);
    }
    v = v/sz;
    //printf("%e \n",v);
    
    return v;
}

//Row-wise variance of the elements of a matrix
void var(double* var_vec, double* M, int nrg, int ncol, double dx)
{
    //double *var_vec;
    double *v1;
    double m;
    //var_vec = (double*) malloc(nrg*sizeof(double));
    v1 = (double*) malloc(ncol*sizeof(double));
    for(int j=0;j<nrg;j++)
    {
    	//print_matr(M,nrg,ncol);
        extract_row(v1, M, j, 0, ncol-1, ncol);
        //print_vec(v1,ncol);
        //m = mean(v1, ncol);
        var_vec[j] = var(v1, ncol, dx);
    }
    
    return;
}

void var(double* var_vec, double* M, int rgstart, int rgend, int colstart, int colend, int sz, double dx)
{
    //double *var_vec;
    double *v1;
    double m;
    //var_vec = (double*) malloc(nrg*sizeof(double));
    v1 = (double*) malloc((colend - colstart +1)*sizeof(double));
    for(int j=rgstart;j<= rgend;j++)
    {
    	//print_matr(M,nrg,ncol);
        extract_row(v1, M, j, colstart, colend, sz);
        //print_vec(v1,ncol);
        //m = mean(v1, ncol);
        var_vec[j-rgstart] = var(v1, colend-colstart+1, dx);
    }
    
    return;
}

double var_complex(fftw_complex* M, int sz)
{
	double v = 0.0;
	for (int i=0;i<sz;i++)
	{
		v += creal(M[i]*conj(M[i]));
	}
	v /= sz;
	
	return v;


}


void var(double* varx, double* vary, double  &var, double* M, int rgstart, int rgend, int colstart, int colend, int sz)
{
	double Sx = 0.0;
	double S2x = 0.0;
	double Sy = 0.0;
	double S2y = 0.0;
	double S = 0.0;
	double S2 = 0.0;
	double c = 1.0/(colend-colstart+1);
	for(int i=rgstart; i<= rgend; i++)
	{
		for(int j=colstart; j<=colend;j++)
		{
			Sx += M[i*sz + j];
			S2x += M[i*sz+j]*M[i*sz+j];
			Sy += M[j*sz + i];
			S2y += M[j*sz+i]*M[j*sz + i];
			S += M[i*sz + j];
			S2 += M[i*sz + j]*M[i*sz + j];
		}
		varx[i-rgstart] = c*S2x - c*c*Sx*Sx;
		vary[i-rgstart] = c*S2y - c*c*Sy*Sy;
		Sx = 0.0;
		S2x = 0.0;
		Sy = 0.0;
		S2y = 0.0;
	}
	var = (c*c*S2 - c*c*c*c*S*S);
	S = 0.0;
	S2 = 0.0;
	
	
	return;
	
	
}

void rescale(double* M, int sz, double fact)
{
	for (int i=0;i < sz;i++)
	{
		M[i] = M[i]/fact;
		
	}
	
	return;
}

void rescale_matrix(double* M, int sz, double fact)
{
	for (int i=0;i < sz;i++)
	{
		for(int j=0;j<sz;j++)
		{
			//printf("M before = %e \n",M[i*sz+j]);
			M[i*sz + j] *= fact;
			//printf("M after = %e \n",M[i*sz+j]);
		}
	}
	
	return;
}



void rescale(fftw_complex* M, int sz, int fact)
{
	for (int i=0;i < sz;i++)
	{
		M[i] = creal(M[i])/double(fact) + I*cimag(M[i])/double(fact);
		
	}
	
	return;
}

void transpose(double* M1, double* M, int sz)
{
	for(int i=0;i < sz; i++)
	{
		for (int j=0; j< sz; j++)
		{
			M1[i*sz+j] = M[j*sz+i];
		}
	}
	
	return;
		
}


//Returns the sum of two vectors
void vec_sum(double* z, double* v1, double* v2, int sz)
{
    //double* z;
    //z = (double*) malloc(sz*sizeof(double));
    for(int i=0;i<sz;i++)
    {
        z[i] = v1[i] + v2[i];
    }
    
    return;
    
}

//Returns the pointwise product of two vectors
void vec_product(double* z, double* v1, double* v2, int sz)
{
    //double* z;
    //z = (double*) malloc(sz*sizeof(double));
    for(int i=0;i<sz;i++)
    {
        z[i] = v1[i]*v2[i];
    }
    
    return;
    
}

void vec_product(fftw_complex* z, fftw_complex* v1, fftw_complex* v2, int sz)
{
    //double* z;
    //z = (double*) malloc(sz*sizeof(double));
    for(int i=0;i<sz;i++)
    {
        z[i] = v1[i]*v2[i];
    }
    
    return;
    
}


//Divide a vector by a scalar
void vec_divide(double* z,double* v1, double n, int sz)
{
    //double* z;
    //z = (double*) malloc(sz*sizeof(double));
    for(int i=0;i<sz;i++)
    {
        z[i] = v1[i]/n;
    }
    
    return;
    
}


//Creates the vector of the modes q
void create_q(double* q, int L, double dx)
{
	//int LL = int(L/dx);
	int Lc = L/2+1;
	
	for(int i=0;i<Lc;i++)
	{
		q[i] = 2*pi*i/(double)L/dx;
	}
	for(int i=Lc;i<L;i++)
	{
		q[i] = 2*pi*(i-L)/(double)L/dx;
	}
	
	
	return;
}

void parameters_ps(double* sigma, double* esigma, double* V, double* qx, double* qy, double nux, double nuy, double kx, double ky, double kxy, double dx, double dt, double D, int L)
{
	int L2 = L*L;
	int Lc = L/2+1;
	
	for (int rg=0; rg<L; rg++)
    {
		for(int col=0; col < Lc;col++)
		{
			//q = sqrt(qx[col]*qx[col] + qy[rg]*qy[rg]);
			sigma[rg*Lc +col] = -( nux*qx[col]*qx[col] + nuy*qy[rg]*qy[rg] + kx*qx[col]*qx[col]*qx[col]*qx[col] + ky*qy[rg]*qy[rg]*qy[rg]*qy[rg] + kxy*qx[col]*qx[col]*qy[rg]*qy[rg] ); 
			esigma[rg*Lc + col] = exp(sigma[rg*Lc+col]*dt);
			V[rg*Lc+col] = sqrt( (D/dx/dx)*(esigma[rg*Lc+col]*esigma[rg*Lc+col] - 1.)/(sigma[rg*Lc+col]+1e-16) );
		}
	}
	
	V[0]=sqrt(2*D*dt/dx/dx); 
    V[L/2]=V[0]; V[L/2*Lc]=V[0]; V[L/2*Lc+L/2]=V[0];
	
	return;
}

void parameters_kdV1d_ps(fftw_complex* sigma, fftw_complex* esigma, double* qx, double omega, double dt, int L)
{
	int L2 = L*L;
	int Lc = L/2+1;
	double x = 0.0;
	double y = 0.0;
	
	for(int col=0; col < Lc;col++)
	{
			//q = sqrt(qx[col]*qx[col] + qy[rg]*qy[rg]);
			y = -omega*qx[col]*qx[col]*qx[col];
			sigma[col] = I*y; 
			esigma[col] = cos(y*dt) + I*sin(y*dt);
			
	}

	
	return;
}

void fourier_transform(FILE *fin, FILE *fout, int Lc, int L, int ntimes)
{
	fftw_complex* Xq = (fftw_complex*) fftw_malloc(Lc*sizeof(fftw_complex));
	double* Xr = (double*) fftw_malloc(L*sizeof(double));
	fftw_plan Xq2Xr = fftw_plan_dft_c2r_1d(L,Xq,Xr,FFTW_ESTIMATE);
		
	float x,y;
	x= 0;
	y=0;
	
	for (int j=0;j<ntimes;j++)
	{
		
		for (int i=0;i<Lc;i++)
		{
			fscanf(fin,"%f+%fi ",&x,&y);
			Xq[i] = x + I*y;
		} 
		
		fftw_execute(Xq2Xr);
		if (j==0)
		{
			rescale(Xr,L,L);
		}
		print_vec2file(Xr, 0, L-1, fout);
		zeros(Xq, Lc);
		zeros(Xr, L);
	}	
	
		
	return;
}
void calculate_psd1d_ps(double* psd, fftw_complex *hqaux, int N, double dx)
{
	int Nc = N/2+1;
	double C = double(1./N/N);
	
	for(int col=0; col<Nc;col++)
	{
		//saux[rg*Lc + col] = cabs(hqaux[rg*Lc+col])*cabs(hqaux[rg*Lc+col]);//C*creal(hqaux[rg*Lc + col]*conj(hqaux[rg*Lc + col]));
		psd[col] = C*creal(hqaux[col]*conj(hqaux[col]));
	}
	return;
}

void calculate_psd2d_ps(double* saux, fftw_complex *hqaux, int L, int Lc, double dx)
{
	//double C = (dx*dx*dx/((double) L*(double) L));
	//double C = double(1./L/L);
	//printf("%e \n",C);
	for(int rg=0; rg<L;rg++)
	{
		for(int col=0; col<Lc;col++)
		{
			saux[rg*Lc + col] = cabs(hqaux[rg*Lc+col])*cabs(hqaux[rg*Lc+col]);//C*creal(hqaux[rg*Lc + col]*conj(hqaux[rg*Lc + col]));
			//saux[rg*Lc + col] = C*creal(hqaux[rg*Lc + col]*conj(hqaux[rg*Lc + col]));
		}
	}
	
	return;
}

void calculate_psd2d(double* saux, fftw_complex *hqaux, int L, double dx)
{
	int LL = int(L/dx);
	int ncol = ((int)LL/2+1);
	double C = (dx*dx*dx/((double) L*(double) L));
	
	for(int rg=0; rg<LL;rg++)
	{
		for(int col=0; col<ncol;col++)
		{
			saux[rg*ncol + col] = C*creal(hqaux[rg*ncol + col]*conj(hqaux[rg*ncol + col]));
		}
	}
	
	
	return;
}

void calculate_psd(double* saux, fftw_complex* hq, int L, double dx, int nvolte)
{
	int LL = int(L/dx);
	double C = (dx*dx/(double) L);
	for(int rg=0;rg<nvolte;rg++)
	{
		//for(int col=1;col<floor((double)(LL/2))+1;col++)
		for(int col=0;col<floor((double)(LL/2));col++)
		{
			//saux[rg*L + col] = (1/(double)L)*hq[rg*L + col]*hq[rg*L + L-col-1];
			//saux[rg*L + col] = (1/(double)L)*cabs(hq[rg*L + col] )*cabs(hq[rg*L + col]);
			//saux[rg*L + col] = (1/(double)L)*hq[rg*L + col]*conj(hq[rg*L + col]);
			
			
			saux[rg*LL + col] = C*creal(hq[rg*LL + col]*conj(hq[rg*LL + col]));
			//saux[rg*L + col] = creal(hq[rg*L + col]*conj(hq[rg*L + col]));
			
			
			//printf("%e + i %e  \n", creal(hq[rg*L + col]), cimag(hq[rg*L + col]));
			//printf("%e + i %e  \n", creal(conj(hq[rg*L + col])), cimag(conj(hq[rg*L + col])));
			//printf("Prod = %e + i %e  \n", creal( hq[rg*L + col]*conj(hq[rg*L +col]) ), cimag(conj(hq[rg*L + col])*hq[rg*L + col]) );
			//saux[rg*L + col] = cabs(hq[rg*L + col] );
			
			//printf("(%e + i%e)*(%e +i%e) \n",creal(hq[rg*L + col]),cimag(hq[rg*L + col]), creal(hq[rg*L + L-col]),cimag(hq[rg*L + L-col]));
			//printf("%e \n",saux[rg*L + col]);
			//printf("%e \n",cimag((1/(double)L)*hq[rg*L + col]*hq[rg*L + L-col]));
		}
	}
	
	return;
}








