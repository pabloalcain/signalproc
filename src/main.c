/* 
 * main routine to interface with the VMD specden plugin 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "specden.h"

int main(int argc, char *argv[]) 
{
  FILE *fp;
  double *input, *output;
  int ndat, nn, specr, normtype;
  double maxfreq, deltat, temp;
  double x, y, z;
  double avg[3];
  int i;
  int n;

  if (argc != 2) {
    fprintf(stderr, "usage: %s fname\n", argv[0]);
    return -1;
  }
  
  fp = fopen(argv[1], "r");
  if (!fp) {
    fprintf(stderr, "Could not read file %s\n", argv[1]);
    return -1;
  }

  n = fscanf(fp, "%d, %d, %d, %lf\n", &ndat, &nn, &specr, &maxfreq);
  if (n != 4) 
    fprintf(stderr, "Wrong format in header of file %s", argv[1]);
  /* defaults */
  normtype = NORM_HARMONIC;
  temp  = 300;
  deltat = 0.01;

  /* size of output data */
  nn = (int) ((double)ndat)*maxfreq/219474.0*deltat/(2.0*M_PI);
  
  if (nn < 1) {
    fprintf(stderr, "Check input file, no data generated!\n");
  }
  
  input = (double *) malloc(3 * ndat * sizeof(double));
  output = (double *) malloc(2 * nn * sizeof(double));

  avg[0]=avg[1]=avg[2]=0.0;
  for (i=0; i<ndat; i++) {
    n = fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);
    if (n != 3) 
      fprintf(stderr, "Wrong format in file %s:%d", argv[1], i+1);
    input[3*i+0] = x;
    input[3*i+1] = y;
    input[3*i+2] = z;
    avg[0] += x;
    avg[1] += y;
    avg[2] += z;
  }
  fclose(fp);

  avg[0] /= (double) ndat;
  avg[1] /= (double) ndat;
  avg[2] /= (double) ndat;
  for (i=0; i<ndat; ++i) {
    input[3*i+0] -= avg[0];
    input[3*i+1] -= avg[1];
    input[3*i+2] -= avg[2];
  }
  
  nn=calc_specden(ndat, input, output, normtype, specr, maxfreq, deltat, temp);

  fp = fopen("output.dat", "w");
  if (!fp) fprintf(stderr, "Could not open file output.dat");
  for (i = 1; i < nn; i++) {
    fprintf(fp, "%g %g\n", output[2*i], output[2*i+1]);
  }
  fclose(fp);
  
  return 0;
}
