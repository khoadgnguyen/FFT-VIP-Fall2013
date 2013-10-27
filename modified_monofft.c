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
#include "fft_sender.h"
#include "fftw3.h"

void test02 (int N, double (*input)[2], double (*output)[2]);


void genfft(float* fbuffer, int N)
{
    printf("IN: monofft:genFFT(): sizeof(fbuffer)=%d\n", sizeof(fbuffer));

    /* Check size of buffer */
    if(sizeof(fbuffer)/sizeof(float) != N)
        err("fbuffer does not contain enough samples");

    int i=0,j=0,k=0,l=0;

    double (*X)[2];                  /* pointer to frequency-domain samples */   
    double x[N][2];             
    double mag, min, max, mag_norm;                 
    X = malloc(2 * N * sizeof(double));  
    while (i<N) {
        x[i][1]=0;
        i++;
    }

    int lSize = sizeof(float)*N; //position of the current stream (It's at the end because of fseek())
    //rewind (f);
	
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

	memcpy(mywave, fbuffer, lSize);	

    int mask=65535;
    int  leftVa;
	while (i<N) {    
	//while (i<((lSize-40)/4)) {

        leftVa=(mywave->data[i+1]<<16) | (mywave->data[i] & mask);
        mywave->data[j]=leftVa;

        /* LEFT CHANNEL DATA GOES HERE */
        x[k][0]=leftVa;

        j++;
        k++;
        //~ printf("k is %d and N is %d\n", k, N);
        if (k==N) {
	    // fft(N,x,X);

		//testing fftw3 here///
            test02(N, x, X);
	
		
            min = 0;
            max = sqrt((X[0][0]*X[0][0])+(X[0][1]*X[0][1]));
            for(l=0; l<N; l++) {
                mag=sqrt((X[l][0]*X[l][0])+(X[l][1]*X[l][1]));
                if(mag > max) {
                    max = mag;
                }
                if(mag < min) {
                    min = mag;
                }
            }


            for(l=0; l<N; l++) {
                mag=sqrt((X[l][0]*X[l][0])+(X[l][1]*X[l][1]));

                mag_norm = (mag - min)/(max - min);
                fbuffer[l] = mag_norm;

            }
            k=0;
        }
        i=i+2;
    }


    free(X);
    return;
}


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
  int n = 100;
  int nc;
  fftw_complex *out;
  fftw_plan plan_backward;
  fftw_plan plan_forward;
  unsigned int seed = 123456789;

/*
  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Demonstrate FFTW3 on a single vector of real data.\n" );
  printf ( "\n" );
  printf ( "  Transform data to FFT coefficients.\n" );
  printf ( "  Backtransform FFT coefficients to recover data.\n" );
  printf ( "  Compare recovered data to original data.\n" );
*/
/*
  Set up an array to hold the data, and assign the data.
*/
  in = fftw_malloc ( sizeof ( double ) * n );

  srand ( seed );

  for ( i = 0; i < n; i++ )
  {
    in[i] = frand ( );
  }

  /*printf ( "\n" );
  printf ( "  Input Data:\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %12f\n", i, in[i] );
  }*/



/*
  Set up an array to hold the transformed data,
  get a "plan", and execute the plan to transform the IN data to
  the OUT FFT coefficients.
*/
  nc = ( n / 2 ) + 1;

  out = fftw_malloc ( sizeof ( fftw_complex ) * nc );

  plan_forward = fftw_plan_dft_r2c_1d ( n, in, out, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );

  printf ( "\n" );
  printf ( "  Output FFT Coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < nc; i++ )
  {
    printf ( "  %4d  %12f  %12f\n", i, out[i][0], out[i][1] );
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
  fftw_destroy_plan ( plan_backward );

  fftw_free ( in );
  fftw_free ( in2 );
  fftw_free ( out );

  return;
}
/******************************************************************************/



