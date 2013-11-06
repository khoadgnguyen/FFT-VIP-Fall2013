/*

    Convert a Stereo file into a Mono wave file by concatenating the left channel data together.
    Passes determined values from wav file into FFT function.
    Normalizes values based on the range, outputted values are from 0-1.
    Writes output to text file for use with the spectrogram display code.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>    
#include "fft.h"
#include <fftw3.h>

/********************** Global Variable ***************************************/
float hamWindow_Multiplier;

/*****************************new added***************************************/
void test02 ( int N, double (*x)[2], double (*X)[2])

/******************************************************************************/
/*
  Purpose:

    TEST02: apply FFT to real 1D data.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 October 2005

  Author:

    John Burkardt
*/
{
  int i;
  double *in;
  double *in2;
  //int n = 100;
  int nc;
  fftw_complex *out;
  fftw_plan plan_backward;
  fftw_plan plan_forward;
  unsigned int seed = 123456789;

/*
  Set up an array to hold the data, and assign the data.
*/
  in = fftw_malloc ( sizeof ( double ) * N );

  
  for ( i = 0; i < N; i++ )
  {
    in[i] = x[i][0];
	
    hamWindow_Multiplier = 0.5*(1-cos(2*3.1415*i/N));	// Multiplying with hamming window
    
    in[i] = hamWindow_Multiplier*in[i];
  }


/*
  Set up an array to hold the transformed data,
  get a "plan", and execute the plan to transform the IN data to
  the OUT FFT coefficients.
*/
  nc = ( N / 2 ) + 1;

  out = fftw_malloc ( sizeof ( fftw_complex ) * nc );

  plan_forward = fftw_plan_dft_r2c_1d ( N, in, out, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );

  
  /*printf ( "\n" );
  printf ( "  Output FFT Coefficients:\n" );
  printf ( "\n" );
 */


  for ( i = 0; i < nc; i++ )	// Saving fftw3 output to X
  {
    X[i][0] = out[i][0];	
    X[i][1] = out[i][1];
  }

/**********************Don't need to recover input*******************

  Set up an arrray to hold the backtransformed data IN2,
  get a "plan", and execute the plan to backtransform the OUT
  FFT coefficients to IN2.

  in2 = fftw_malloc ( sizeof ( double ) * n );

  plan_backward = fftw_plan_dft_c2r_1d ( n, out, in2, FFTW_ESTIMATE );

  fftw_execute ( plan_backward );

  printf ( "\n" );
  printf ( "  Recovered input data divided by N:\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %12f\n", i, in2[i] / ( double ) ( n ) );
  }
*/


/*
  Release the memory associated with the plans.
*/
  fftw_destroy_plan ( plan_forward );
  //fftw_destroy_plan ( plan_backward );

  fftw_free ( in );
  //fftw_free ( in2 );
  fftw_free ( out );

  return;
}
/******************************************************************************/

void genfft(float* fbuffer, int N)
{
    int i=0,j,k,l=0;
    int wavbuffer[N];
    double (*X)[2];                  /* pointer to frequency-domain samples */   /* double */
    double x[N][2];             /* double */
    double mag, min, max, mag_norm;                 /* double */
    X = malloc(2 * N * sizeof(double));  /* double */
    while (i<N) {
        x[i][1]=0;
        i++;
    }

    mkfifo("wavePipe", S_IRWXO);
    system("arecord -d1 -r44100  -D plughw:CARD=Snowflake,DEV=0  wavePipe");

    int lSize;
    FILE * f = fopen("wavePipe", "r"); //opening the 2 channels wave file

    if(!f) printf("Error reading from wavePipe");

    //just to see what we're dealing with
    fseek (f , 0 , SEEK_END);
    lSize = ftell (f);
    rewind (f);

    typedef struct wave {
        int chunkID;
        int chunkSize;
        int format;
        int subChunk1ID;
        int subChunk1size;
        int audio_num_ch;
        int sampleRate;
        int byteRate;
        int bckAli_bitPs;
        int subChunk2ID;
        int subChunk2size;
        int  data[(lSize-40)/4];
    } wave;
    wave * mywave;
    mywave=(wave *) malloc( lSize);

    fread(mywave,1,lSize,f);

    i=0,j=0;
    k=0;
    int mask=65535;
    int  leftVa;
    while (i<((lSize-40)/4)) {

        leftVa=(mywave->data[i+1]<<16) | (mywave->data[i] & mask);
        mywave->data[j]=leftVa;

        x[k][0]=leftVa;
        wavbuffer[k] = leftVa;
        j++;
        k++;

        if (k==N) {
            test02(N, x, X);
            min = 0;
            max = sqrt((X[0][0]*X[0][0])+(X[0][1]*X[0][1]));
            for(l=0; l<((N/2)+1); l++) {
                mag=sqrt((X[l][0]*X[l][0])+(X[l][1]*X[l][1]));
                if(mag > max) {
                    max = mag;
                }
                if(mag < min) {
                    min = mag;
                }
            }


            for(l=0; l<((N/2)+1); l++) {
                mag=sqrt((X[l][0]*X[l][0])+(X[l][1]*X[l][1]));

                mag_norm = (mag - min)/(max - min);
                fbuffer[l] = mag_norm;

            }
            k=0;
        }
        i=i+2;
    }


    
    //for(i = 0; i < N; i++)
    //    printf("wavbuffer[%d] is %12d, fbuffer[%d] is %12f\n", i, wavbuffer[i], i, fbuffer[i]);


   // fwrite(mywave,1,(lSize/2)+22,f1);
    fclose (f);
    // fclose (f2);
    free(X);
    return;
}



