#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
//#include <iostream>
//#include <malloc.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


//#include "gauss.h"
#include "/home/edoardo/Documents/UC3M/Dottorato/Equations/equations.h"

const gsl_rng_type * T;
gsl_rng * r;

    



//Reads parameters from file
void read_params(double *nu, double *lambda, double *D, double *dx, double *dt, int *N, int *tmax, double *tmin, int *nstat, int *nvolte, int *ntimes_psd, FILE* fparam)
{
    char s[20];
    fscanf(fparam,"%s",s); fscanf(fparam,"%le",nu);
    fscanf(fparam,"%s",s); fscanf(fparam,"%le",lambda);
	fscanf(fparam,"%s",s); fscanf(fparam,"%le",D);
    fscanf(fparam,"%s",s); fscanf(fparam,"%le",dx);
    fscanf(fparam,"%s",s); fscanf(fparam,"%le",dt);
    fscanf(fparam,"%s",s); fscanf(fparam,"%d",N);
    fscanf(fparam,"%s",s); fscanf(fparam,"%d",tmax);
    fscanf(fparam,"%s",s); fscanf(fparam,"%le",tmin);
    fscanf(fparam,"%s",s); fscanf(fparam,"%d",nstat);
    fscanf(fparam,"%s",s); fscanf(fparam,"%d",nvolte);
    fscanf(fparam,"%s",s); fscanf(fparam,"%d",ntimes_psd);
    //printf("%le \n",*nu);
    return;
}

void gen_ns(double *noise, int sz, gsl_rng * r)
{
	for (int j=0;j<sz;j++)
	{
		//noise[j] = randn_notrig();
		noise[j] = gsl_ran_gaussian_ziggurat(r, 1.0);
		//noise[j] = gsl_ran_flat(r,-0.5,0.5);
		//noise[j] = gsl_ran_gaussian_ziggurat(r, 1.0);
		
	}
	
	return;
}

//Initial conditions
/*
void initial_conditions(double* V, int N, double dx, double lambda, double omega, double v)
{
	double C = 3*omega*v/lambda;
	double C1 = sqrt(v)/2;
	for(int i=0;i<N;i++)
	{
		V[i] = C/cosh(C1*(i-N/10)*dx)/cosh(C1*(i-N/10)*dx);
		//V[i] = C/cosh(C1*(i)*dx)/cosh(C1*(i)*dx);
	}
	//V[0] = V[N-1];
	
	return;
}
*/


//Padding
void padd(fftw_complex* aq, fftw_complex* hq, fftw_complex* bq, double* qx, int M, int Mc, int N, int Nc)
{
	//
	zeros(aq, Mc);
	zeros(bq, Mc);
	
	for(int i=0;i<Nc;i++)
	{
		aq[i] = hq[i];
		bq[i] = I*(qx[i]*hq[i]);
	}
	

		
	return;
}

/*
void convolution(fftw_complex* aq, fftw_complex* hq, fftw_complex* bq, double * qx,int N, int Nc)
{
	zeros(aq, Nc);
	for(int i=0;i<Nc;i++)
	{
		bq[i] = I*qx[i]*hq[i];
	}
	for(int i=0;i<Nc;i++)
	{
		for(int n=0;n<Nc;n++)
		{
			aq[n] += bq[i]*hq[n-i];
		}
	}

	
	
	return;
}
*/

//Main Loop Integration
void integrate_ps(double* h, fftw_complex* hq,double* h_aux, fftw_complex* hq_aux,  double* V, double* ar, fftw_complex* aq, double* br, fftw_complex* bq, double* b,double* sigma, double* esigma, double* qx, double* t, double dt, double* delta_t, int nstat, int ntimes_psd, int N, double dx, double nu, double lambda, FILE *fout_h, FILE *fout_w)
{
    int k, rg, col, nnr, nnl;
   int LL = int(N/dx);
    int LLaug = LL+4;
    
    int L2 = N*N;
    int Nc = int(N/2 + 1);
    int M = int(3*N/2);
    int Mc = int(M/2+1);
    
    double fact = double(N)/double(M);
    double w = 0.0;
    
    
    int npsd = int(nstat/ntimes_psd);
    int cnt = 0;
    int x = 0;
    char str[10];
    FILE *fout_s;
    
  
    register int i;
    
    
    
    printf("N= %d \n", N);
    
    double* noise = (double*) fftw_malloc(N*sizeof(double));
    fftw_complex* noiseq = (fftw_complex*) fftw_malloc(Nc*sizeof(fftw_complex));
    
    double* psd = (double*) fftw_malloc(Nc*sizeof(double));
   
   
    fftw_plan hq2h, h2hq, ar2aq, aq2ar, bq2br, br2bq, hq_aux2h, h2hq_aux, hq_aux2h_aux, noise2noiseq;

    h2hq = fftw_plan_dft_r2c_1d(N,h,hq,FFTW_ESTIMATE);

    hq2h = fftw_plan_dft_c2r_1d(N,hq,h,FFTW_ESTIMATE);
    
    
    hq_aux2h = fftw_plan_dft_c2r_1d(N,hq_aux,h,FFTW_ESTIMATE);
    hq_aux2h_aux = fftw_plan_dft_c2r_1d(N,hq_aux,h_aux,FFTW_ESTIMATE);
    h2hq_aux = fftw_plan_dft_r2c_1d(N,h,hq_aux,FFTW_ESTIMATE);
	
    ar2aq = fftw_plan_dft_r2c_1d(M, ar, aq, FFTW_ESTIMATE);

    aq2ar = fftw_plan_dft_c2r_1d(M, aq, ar, FFTW_ESTIMATE);
    
    bq2br = fftw_plan_dft_c2r_1d(M, bq, br, FFTW_ESTIMATE);
    br2bq = fftw_plan_dft_r2c_1d(M, br, bq, FFTW_ESTIMATE);
    
    noise2noiseq = fftw_plan_dft_r2c_1d(N, noise, noiseq, FFTW_ESTIMATE);
    
    
  
    
   //Initial condition
   
   /*
   FILE *f_init;
   f_init = fopen("initial.dat","w");
   printf("Initial conditions...\n");
   initial_conditions(h, N, dx, lambda, omega, v);
   print_vec2file(h,0,N-1, f_init);
   fclose(f_init);
   fftw_execute(V2Vq);
   //rescale(hq, Nc, N);
   //rescale(hq, Nc, sqrt(N));
   printf("Initial conditions...done\n");
   * 
   */
   zeros(h, N);
   zeros(hq, Nc);
   zeros(ar, N);
   zeros(br, N);
   zeros(aq, Nc);
   zeros(bq, Nc);
   zeros(h_aux, N);
   zeros(hq_aux, Nc);
   zeros(psd, Nc);
   
   
   double C1 = lambda*fact*dt;
 
   //fout_V = fopen("h.dat","w");
   for (k=0; k < nstat; k++)
    {
		printf("%e \n",t[k]);
		fflush(stdout);
		//for (int z=0; z< npsd; z++)
		//{
	       	for(i=0;i<delta_t[k];i++)
	       	{
				gen_ns(noise, N, r);
				fftw_execute(noise2noiseq);
				
	       		padd(aq, hq, bq, qx, M, Mc, N, Nc);
	       		//printf("Transform..\n");
	       		fftw_execute(aq2ar);
	       		fftw_execute(bq2br);
	       		//printf("Transform..done\n");
				vec_product(ar, ar, br, M);
				//fftw_execute(ar2aq);	       		
				//fftw_execute(V2Vq);
	       		rescale(ar, M, double(N*N));
	       		
	       		fftw_execute(ar2aq);
	       		//printf("Transform..done\n");
	       		
	       		
	       		
	       		
	       		for(col=0; col< Nc; col++)
				{
					
					
					//hq[col] = hq[col] + dt*(sigma[col]*hq[col] + lambda*fact*aq[col]) + V[col]*noiseq[col];
					hq[col] = esigma[col]*(C1*aq[col] + hq[col]) + V[col]*noiseq[col];
					
				
				}
				
				               
	        }
	        
	        cnt++;
			if (cnt == npsd)
			{
				x = sprintf(str,"psd%e.dat",t[k]);
				fout_s = fopen(str,"a+");
		
				
				calculate_psd1d_ps(psd, hq, N, dx);
				//print_mat2file(saux, N, Nc, fout_s);
				print_vec2file(psd, 0, Nc-1, fout_s);
				fclose(fout_s);
				cnt = 0;
			}

	        //printf("copy \n");
	        copy_vec(hq, hq_aux, Nc);
	        
			//printf("transform \n");
	        fftw_execute(hq_aux2h_aux);
	        
	        //printf("rescale \n");
			rescale(h_aux, N, N);
					
	       
	        print_vec2file(h_aux,0,N-1, fout_h);
	        
	        //fftw_execute(V2Vq);
	        
	
	        //copy_vec(hq_aux,hq,Nc);
	        
	        
	        
	      
	        //printf("Anti-Transform..\n");
	        //fftw_execute(V2Vq);
	        //printf("Anti-Transform..done!\n");
	      
	       
	       
	        
	        //Calculate roughness
	        
	        //fftw_execute(hq2h);
	        //rescale(haux2, N, N*N);
	        //var(roughx, roughy, w, haux2, 0, N-1, 0, N-1, N);
	        //wx = mean(roughx, 0, N-1);
	        //wy = mean(roughy, 0, N-1);
	        w = var(h_aux, N, dx);
	        fprintf(fout_w, "%e ",w);
	        fflush(fout_w);
	        
	        
	        
	        //Print roughness to file
	        /*
	        fprintf(fout_wx,"%e %e \n",t[k],wx);
	        fflush(fout_wx);
	        fprintf(fout_wy,"%e %e \n",t[k],wy);
	        fflush(fout_wy);
	        fprintf(fout_w,"%e %e \n",t[k],w);
	        fflush(fout_w);
			fftw_execute(h2hq);
			
			
			
			
			
       */
    }
    
    fprintf(fout_w, "\n ");
	fflush(fout_w);
    /*
    print_mat2file(haux2, 0, N-1, 0, N-1, N, fout_h);
    //Calculate psd
    fout_s = fopen("psd.dat","a+");
	
	calculate_psd2d_ps(saux, hqaux2, N, Nc, dx);

	//Print psd to file

	print_mat2file(saux, N, Nc, fout_s);
	fflush(fout_s);
	fclose(fout_s);
	
  */
	
    
    return;
    
}

void parameters_HK1d_ps(double* V, double* sigma, double* esigma, double* qx, double nu, double D, double dt, double dx, int N)
{
	//int L2 = L*L;
	int Nc = N/2+1;
	double x = 0.0;
	double y = 0.0;
	
	for(int col=0; col < Nc;col++)
	{
			//q = sqrt(qx[col]*qx[col] + qy[rg]*qy[rg]);
			sigma[col] = -nu*qx[col]*qx[col];
			esigma[col] = exp(sigma[col]*dt);
			V[col] = sqrt( (D/dx)*(esigma[col]*esigma[col]-1.)/(sigma[col]+1e-16) );
			
	}
	V[0] = sqrt(2*D*dt/dx);
	V[Nc-1] = V[0];

	
	return;
}


int main(int narg,char **args)
{
	
    srand(time(0));
    
    gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	
	gsl_rng_set (r, time(0));
    
    
    //Height of the surface; it's a matrix of nvolte rows and N columns; haux are auxiliar vectors
   // double *h, *Vaux3;//*haux1, *haux2, *haux3;
    double *s, *saux, *saux2;
    double *qx, *qy;
    double *sigma, *esigma, *V;//, *W;
    double *noise;
    //fftw_complex *hq;
    fftw_complex *Vqaux;
    fftw_complex *noiseq;
    //fftw_plan p, pnoise, Vq2V;
    
    //Parameters
    int N, Laug, tmax, nstat, nvolte, ntimes_psd;
    double lambda, nu, dx, dt, tmin, D;
    int LL, LLaug;
    int L2, Nc, M, Mc;
    FILE* fparam;
    //Times
    double *times;
   // double tmin;
    double *delta_t; //Interval beetween two consecutive times, divided by dt and approximated to integer. It gives me how many integration loops are there beetween one point and the next one
    
    
    //Observables
    double *w; //Roughness
    double *w1;
    
    //Output files
    FILE *fout_s, *fout_h, *fout_Vreal, *fout_w;
    //fout_wx, *fout_wy,
    //Reading parameters from file
    fparam = fopen(args[1],"r");
    read_params(&nu, &lambda, &D, &dx, &dt, &N, &tmax, &tmin,&nstat, &nvolte, &ntimes_psd, fparam);
    fclose(fparam);
    printf("%d \n",N);
    fflush(stdout);
    
    LL= int(N/dx);
    L2 = N*N;
    Nc = N/2+1;
    M = int(3*N/2);
    Mc = int(M/2+1);
   
   //Initializing times
    times = (double*) malloc(nstat*sizeof(double));
    delta_t = (double*) malloc(nstat*sizeof(double));
    logtimes(times, tmin, tmax, dt, delta_t, nstat);
    
    FILE* f_times = fopen("times.dat","w");
    print_vec2file(times, 0, nstat-1, f_times);
    fclose(f_times);
    
    
    //Initializing the surface
   
    qx = (double*) fftw_malloc(N*sizeof(double));
   // qy = (double*) fftw_malloc(N*sizeof(double));
    
    sigma = (double*) fftw_malloc(N*sizeof(double));
    esigma = (double*) fftw_malloc(N*sizeof(double));
    V = (double*) fftw_malloc(N*sizeof(double));
   // W = (double*) fftw_malloc(N*Nc*sizeof(double));

   
    
    create_q(qx,N,dx);
    FILE *f_q = fopen("qx.dat","w");
    print_vec2file(qx,0,N-1,f_q);
    fclose(f_q);
    //create_q(qy,N,dx);
     
     printf("Initialize variables...\n");
     double* h = (double*) fftw_malloc(N*sizeof(double));
    fftw_complex* hq = (fftw_complex*) fftw_malloc(Nc*sizeof(fftw_complex));
    double* h_aux = (double*) fftw_malloc(N*sizeof(double));
    fftw_complex* hq_aux = (fftw_complex*) fftw_malloc(Nc*sizeof(fftw_complex));
    
    
    double* ar = (double*) fftw_malloc(M*sizeof(double));
    fftw_complex* aq = (fftw_complex*) fftw_malloc(Mc*sizeof(fftw_complex));
    double* b = (double*) fftw_malloc(M*sizeof(double)); 
    double* br = (double*) fftw_malloc(M*sizeof(double));
    fftw_complex* bq = (fftw_complex*) fftw_malloc(Mc*sizeof(fftw_complex));
    printf("Initialize variables...done\n");
   
   

    parameters_HK1d_ps(V, sigma, esigma, qx, nu, D, dt, dx, N);

    
    
    
    fout_h = fopen("h.dat","w");
    //fout_s = fopen("psd.dat","w");
    //fout_wx = fopen("roughx.dat","w");
    //fout_wy = fopen("roughy.dat","w");
   
    fout_w = fopen("rough.dat","w");
    for(int k=0; k<nvolte; k++)
    {
     	
     	//Main Integration
		printf("Integration number %d begins: \n", k+1);
		
		printf("Using PS integration: \n");
		fflush(stdout);
		
		integrate_ps(h, hq, h_aux, hq_aux, V, ar, aq, br, bq, b, sigma, esigma, qx, times, dt, delta_t, nstat, ntimes_psd, N, dx, nu, lambda, fout_h, fout_w);
		//integrate_ps(h, qx, qy, esigma, W, times, delta_t, nstat, ntimes_psd, N, dx, fout_V, fout_s, fout_wx, fout_wy, fout_w, r);

		printf("Integration number %d ended. \n", k+1);
		fflush(stdout);
         
       
		
		
    }
    fclose(fout_h);
    //fout_V = fopen("h.dat","r");
    //fout_Vreal = fopen("Vreal.dat","w");
    //fourier_transform(fout_V, fout_Vreal, Nc, N, ntimes_psd);
     
    // fclose(fout_V);
     //fclose(fout_s);
     //fclose(fout_wx);
     //fclose(fout_wy);
     fclose(fout_w);
     
    
   
  
    
    
   return 0;
    
}
