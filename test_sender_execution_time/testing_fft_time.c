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
# include <time.h>


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

    //mkfifo("wavePipe", S_IRWXO);
    //system("arecord -d1 -r44100  -D plughw:CARD=WEBCAM,DEV=0  wavePipe");

    int lSize;
    FILE * f = fopen("test.txt", "r"); //opening the 2 channels wave file
	
    if(!f) printf("Error reading from wavePipe");

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
        //~ printf("k is %d and N is %d\n", k, N);
        if (k==N) {
            fft(N, x, X);
	    
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


    
    //for(i = 0; i < N; i++)
	//printf("test: %f , %f \n", x[i][0],X[i][1]);
    //    printf("wavbuffer[%d] is %12d, fbuffer[%d] is %12f\n", i, wavbuffer[i], i, fbuffer[i]);

    fclose (f);
    free(X);
   
    return;
}


