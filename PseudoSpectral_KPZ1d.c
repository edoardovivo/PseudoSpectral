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
#include "equations.h"

const gsl_rng_type * T;
gsl_rng * r;

struct params{
	int N; 							//Number of lattice points
	double nu, lambda, D; 			//The parameters of the KPZ equation
	double dx, dt; 					//Space and time discretization
	double tmax; 					//Final integration time

	int nTimesOut; 					//For each run of the integration, the number of times the
  									//observables are computed and stored.

	int nPSD; 						//For each run of the integration, the number of times the
				  	  	  	  	  	  	//power spectral density is computed and stored.


	int nStat; 						//Number of runs of the integration

};

struct FFTplans{

};

struct Interface{


};

int main(int narg,char **args)
{


	//Initialize random number generator


	//Initialize parameters and auxiliary variables structures

	//Read parameters

	//Initialize FFT plans


	//Initialize vectors of times, momenta, and the surface in momentum space


	//Compute auxiliary variables (sigma, exp(sigma), etc..)


	//Integrate the equation and output the results


	return 0;
}
