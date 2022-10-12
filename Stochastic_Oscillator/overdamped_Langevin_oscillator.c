/***************

Numerical integration of the overdamped Langevin equation for a stochastic oscillator.


"Hands-on" demonstraton for Nordita School in Stockholm 
by C. Micheletti, michelet@sissa.it

version 1.1,  01/03/2012



***************/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
float gasdev(long *idum); /* Numerical recipes routine: Gaussian random numbers */



main(){


  double position, dt, time, gamma, KbT, eta;
  double new_position, mean_position, mean_square_position;
  double k;
  int i;
  long int n_measurements;
  long idum;
  FILE *fp2;

  idum=-1;  /* seed to initialise random number generator */
  KbT=1.0;
  gamma=1.0;
  k=10.0;
  dt = 0.01*gamma/k; /* integration time step, much smaller than gamma/k */

  printf("\n\n\nKbT:   %lf\ngamma: %lf\nk:     %lf\ntau=gamma/k: %lf\nintegration_dt: %lf",KbT,gamma,k,gamma/k,dt);

  mean_position=0;
  mean_square_position=0;
  n_measurements=0;




  fp2 = fopen("position_overdamped_oscillator.dat","w");
  

  position=0.0;
  
  for(i=0; i <= 1000000; i++){ /* advance time step */
    time = i*dt;
    
    /* find velocity and position at next time step */
      
    eta = gasdev(&idum); /* random number. zero mean and unit variance */
    new_position = position -k*position*dt/gamma + sqrt(2*KbT*gamma) * eta * sqrt(dt)/gamma;
    
    position = new_position; 
    
    fprintf(fp2,"%lf %lf\n",time,position); /* write our position */
    
    /* Cumulate the position and square position for later averaging*/
    
    mean_position += position;
    mean_square_position+= position*position; 
    n_measurements++;
  }

  fclose(fp2);

  mean_square_position=mean_square_position/n_measurements;  
  mean_position=mean_position/n_measurements;
  printf("\n\nMean position:         %lf    Expected: %lf\n", mean_position,0.0);  
  printf("Mean square position:  %lf    Expected: %lf\n", mean_square_position, KbT/k);
 
}






#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran2 (long int *idum)
{
  	int j;
  	long k;
  	static long idum2 = 123456789;
  	static long iy = 0;
  	static long iv[NTAB];
  	float temp;

  	if (*idum <= 0)
    {
      	if (-(*idum) < 1)
			*idum = 1;
      	else
			*idum = -(*idum);
      	idum2 = (*idum);
      	for (j = NTAB + 7; j >= 0; j--)
		{
	  		k = (*idum) / IQ1;
	  		*idum = IA1 * (*idum - k * IQ1) - k * IR1;
	  		if (*idum < 0)
	    		*idum += IM1;
	  		if (j < NTAB)
	    		iv[j] = *idum;
		}
      	iy = iv[0];
    }
  	k = (*idum) / IQ1;
  	*idum = IA1 * (*idum - k * IQ1) - k * IR1;
  	if (*idum < 0)
    	*idum += IM1;

  	k = idum2 / IQ2;
  	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  	if (idum2 < 0)
    	idum2 += IM2;
  	j = iy / NDIV;
  	iy = iv[j] - idum2;
  	iv[j] = *idum;
  	if (iy < 1)
    	iy += IMM1;
  	if ((temp = (float)AM * iy) > RNMX)
    	return RNMX;
  	else
    	return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


float gasdev(long *idum)
{
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

