/************************************************************************
 * This program takes time-data and calculates the
 * powerspectrum/fourier transform of the autocorrelation function.
 * This is a re-write of the fourier.x code written by Volker Kleinschmidt
 * and Harald Forbert as a tcl plugin for VMD by Axel Kohlmeyer.
 * (c) 2002-2005 Harald Forbert, Volker Kleinschmidt (c) 2002-2008 Axel Kohlmeyer.
 *
 * usage: calc_specden(<ndat>,<input>,<output>,<deltat>,<maxfreq>,<temp>,<specr>);
 * <ndat>    number of data sets.
 * <input>   time series data.
 * <output>  power spectrum.
 * <normtype> normalization correction type (fourier, classic, kubo, harmonic, schofield)
 * <deltat>  time difference between data sets (in atomic units).
 * <maxfreq> max fequency (in wavenumbers).
 * <temp>    temperature (in kelvin)
 * <specr>   resolution of spectrum (1 gives maximal resolution and noise).
 *
 * the various corrections are:
 * fourier:    is the plain power spectrum of the input data (normalized to
 *             unity in the output frequency range.
 * classical:  is the power spectrum with a prefactor of 
 *             \omega ( 1 - \exp(-\beta \hbar \omega) )
 *             corresponding to the classical/Gordon limit.
 * kubo:       is the power spectrum with a prefactor of
 *             \omega \tanh(\beta \hbar \omega/2)
 *             corresponding to the Kubo correction
 * harmonic:   is the power spectrum with a prefactor of
 *             \omega \beta \hbar \omega
 *             corresponding to the high temperature / harmonic limit
 *             NOTE: this is the _recommended_ correction factor.
 * schofield:  is the power spectrum with a prefactor of
 *             \omega ( 1 - \exp(-\beta \hbar \omega) ) *
 *                                              \exp(\beta \hbar \omega /2)
 *             corresponding to Schofield's correction
 *
 * All spectra with their corresponding prefactor are separately normalized
 * in the output range to sum up to unity.
 *
 * Note: the index of refraction of the medium is set to unity.
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "specden_cu.h"

typedef union {double re; double im;} cmplx;

/* helper function:
 *
 * calculate f_sum   = cos_sum**2 + sin_sum**2,
 *      with cos_sum =  sum_j cos(j*w) input_j
 *       and sin_sum =  sum_j sin(j*w) input_j
 *
 * sums start at 0, but indices of input start at 1, e.g.
 * cos_sum = sum_{j=1}^n cos((j-1)*w) input_j
 */
__device__ static void fourier_sum (const int n, const double *input, const double omega, cmplx *result)
{
  int k;
  double  lambda, duk, uk, cs, ss, cf, sf;
  
  /* in order to be able to sum up the input_j in ascending order
   * use above algorithm with inverse data ordering and correct
   * omega -> - omega and the new origin at the end */

  uk = 0.0;
  duk = 0.0;
  if (cos(omega) > 0.0) {
    lambda = -4.0*sin(0.5*omega)*sin(0.5*omega);
    for (k=1; k <= n; ++k) {
      uk  = uk + duk;
      duk = lambda*uk + duk + input[3*k];
    }
  } else { /* cos(omega) <= 0.0_dbl */
    lambda = 4.0*cos(0.5*omega)*cos(0.5*omega);
    for (k=1; k <= n; ++k) {
      uk  = duk - uk;
      duk = lambda*uk - duk + input[3*k];
    }
  }
  cs = duk - 0.5 * lambda * uk;
  ss = uk * sin(omega);

  /* now correct for ordering: */
  cf = cos(omega*(n-1));
  sf = sin(omega*(n-1));

  result->re = cf*cs+sf*ss;      /* cos_sum */
  result->im = sf*cs-cf*ss;      /* sin_sum */

  return;
/*  return cf*cs*cf*cs+sf*ss*sf*ss+sf*cs*sf*cs+cf*ss*cf*ss; */
}

/* main function */
int calc_specden(const int ndat, double *input, double *output, 
                 const int normtype, const int specr, 
                 const double maxfreq, const double deltat, const double temp) 
{
  int    nn, i;
  double wave_fac, bh, dt, t, f, e;

  double *cu_ftrans, *cu_wtrans, *cu_input;
  double *ftrans, *wtrans;
  double norm_fourier, norm_classic, norm_kubo, norm_harmonic, norm_schofield;

  int bytes;

  wave_fac = 219474.0/deltat;
  bh       = 1.05459e-34/1.38066e-23/2.41889e-17/deltat/temp;
  
  if (specr < 1) {
    fprintf(stderr, "\nspecden spectrum resolution factor must be bigger or equal 1.\n");
    return -20;
  }

  /* number of frequencies */
  nn = (int) ((double)ndat)*maxfreq/wave_fac/(2.0*M_PI);
  if (nn+1 > ndat) {
    fprintf(stderr, "Maximum frequency too large\n");
    return -40;
  }
  nn = nn/specr;
  
  cudaMalloc((void **) &cu_ftrans, (nn+2)*sizeof(double));
  if (cu_ftrans == NULL) {
    fprintf(stderr, "Out of memory, while trying to allocate array 'ftrans'.\n");
    return -50;
  }
  cudaMalloc((void **) &cu_wtrans, (nn+2)*sizeof(double));
  if (cu_wtrans == NULL) {
    fprintf(stderr, "Out of memory, while trying to allocate array 'wtrans'.\n");
    return -60;
  }

  dt = 2.0*specr*M_PI/((ndat+1)*specr);
    
  /* alloc and copy array to device */

  bytes = (ndat + 2)*3*sizeof(double);
  cudaMalloc((void **) &cu_input, bytes);
  cudaMemcpy(cu_input, input, bytes, cudaMemcpyHostToDevice);
  
  /* compute */
  int nblocks = ndat/TPB + 1;
  dim3 grid((nn+2), specr/TPB + 1);
  dim3 block(1, TPB);

  window<<<nblocks, TPB>>>(cu_input, ndat);
  compute<<<grid, block>>>(cu_input, ndat, nn, specr, dt, cu_ftrans, cu_wtrans);
  
  ftrans = (double *) malloc((nn+2)*sizeof(double));
  wtrans = (double *) malloc((nn+2)*sizeof(double));

  cudaMemcpy(ftrans, cu_ftrans, (nn+2)*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(wtrans, cu_wtrans, (nn+2)*sizeof(double), cudaMemcpyDeviceToHost);
  /* compute norm */
  norm_fourier=norm_classic=norm_kubo=norm_harmonic=norm_schofield=0.0;
  for (i=0; i<=nn; ++i) {
    t = wtrans[1+i];
    f = ftrans[1+i];
    e = t*(1.0 - exp(-bh*t));
    
    norm_fourier  += f;
    norm_classic  += f*e;
    norm_kubo     += f*e/(1.0+exp(-bh*t));
    norm_harmonic += f*t*t;
    norm_schofield += f*e*exp(0.5*bh*t);
  }
  norm_fourier  = 1.0/norm_fourier;
  norm_classic  = 1.0/norm_classic;
  norm_kubo     = 1.0/norm_kubo;
  norm_harmonic = 1.0/norm_harmonic;
  norm_schofield = 1.0/norm_schofield;

  /* output */
  for (i=0; i<=nn; ++i) {
    t = wtrans[1+i];
    f = ftrans[1+i];
    e = t*(1.0 - exp(-bh*t));

    output[2*i] = wave_fac*t;
    switch (normtype) {
      case NORM_FOURIER:
         output[2*i+1] = norm_fourier*f;
         break;
      case NORM_CLASSIC:
         output[2*i+1] = norm_classic *f*e;
         break;
      case NORM_KUBO:
         output[2*i+1] = norm_kubo*f*e/(1.0+exp(-bh*t));
         break;
      case NORM_HARMONIC:
         output[2*i+1] = norm_harmonic*f*t*t;
         break;
      case NORM_SCHOFIELD:
         output[2*i+1] = norm_schofield*f*e*exp(0.5*bh*t);
         break;
      default:
         fprintf(stderr, "specden: unknown normalization. %d\n", normtype);
         return -200;
    }
  }
  return nn;
}




    
__global__ void window(double* cu_input, int ndat) {
  int idx;

  idx = blockDim.x * blockIdx.x + threadIdx.x;
  double win;
  
  if (idx != 0 && idx < ndat+2) {
    win=((double)(2*idx-ndat-1))/((double)(ndat+1));
    win=1.0-win*win;
    cu_input[3*idx]   *=win;
    cu_input[3*idx+1] *=win;
    cu_input[3*idx+2] *=win;
  }
}

__global__ void compute(double *cu_input, int ndat, int nn, int specr, double dt, double *cu_ftrans, double *cu_wtrans){
  
  /* So far threadIdx.x will always be 0. The usefulness of this will
     depend on how much we want to average the signal */

  __shared__ float c[1024];
  int i; /* This index runs over the nn data */
  int k; /* This one runs over the "averaging" */
  double f, s, e;
  i = blockDim.x * blockIdx.x + threadIdx.x;
  k = blockDim.y * blockIdx.y + threadIdx.y;

  if (i < nn+1) {
    cmplx f1,f2,f3;

    double t = 2.0*((double)(i*specr))*M_PI/((double)(ndat+1));
    
    if (k < specr) {

      /* sum over all three dimensions */
      fourier_sum(ndat,(cu_input+0), t+(double)k*dt, &f1);
      fourier_sum(ndat,(cu_input+1), t+(double)k*dt, &f2);
      fourier_sum(ndat,(cu_input+2), t+(double)k*dt, &f3);
      f = f1.re*f1.re;
      f += f1.im*f1.im;
      f += f2.re*f2.re;
      f += f2.im*f2.im;
      f += f3.re*f3.re;
      f += f3.im*f3.im;
      
      /* input data should have zero mean... */
      if (i+k == 0) f=0.0;
        
      /* apply cubic spline correction for input data */
      s=0.5*(t+k*dt);
        
      if (s>0.1) {
        e=pow(sin(s)/s,4.0);
      } else {
        e=pow(1.0-(s*s)/6.0+(s*s*s*s)/120.0,4.0);
      }
      e = e*3.0/(1.0+2.0*cos(s)*cos(s));
      c[k] = e*e*f;
    }
    
    for (int s=blockDim.y/2; s>0; s>>=1) {
      if (threadIdx.y < s) {
        c[threadIdx.y] += c[threadIdx.y + s];
      }
      __syncthreads();
    }
    cu_wtrans[1+i] = t+0.5*dt*((double)(specr-1));
    cu_ftrans[1+i] = c[0];
  }
}
