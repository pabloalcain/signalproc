/* common definitions for the specden plugin */
#ifndef _SPECDEN_PLUGIN_H
#define _SPECDEN_PLUGIN_H
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TPB 128

enum specden_norm_type {
  NORM_FOURIER, NORM_CLASSIC, NORM_KUBO, NORM_HARMONIC, NORM_SCHOFIELD
};

__global__ void window(double *cu_input, int ndat);
__global__ void compute(double *cu_input, int ndat, int nn, int specr, double dt, double *cu_ftrans, double *cu_wtrans);

#ifdef __cplusplus
extern "C" 
{
#endif

extern int calc_specden(const int ndat, double *input, double *output,
                        const int normtype, const int specr, 
                        const double maxfreq, const double deltat, 
                        const double temp);
  
extern double *calc_sgsmooth(const int ndat, double *input,
                             const int window, const int order);

extern double *calc_sgsderiv(const int ndat, double *input,
                             const int window, const int order, 
                             const double delta);
#ifdef __cplusplus
}
#endif
