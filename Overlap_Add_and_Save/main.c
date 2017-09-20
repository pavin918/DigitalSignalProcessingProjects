#include "dsk6713_aic23.h"
#include "dsk6713_led.h"

#define DSK6713_AIC23_INPUT_MIC 0x0015
#define DSK6713_AIC23_INPUT_LINEIN 0x0011

Uint32 fs = DSK6713_AIC23_FREQ_16KHZ;                           // 1
Uint16 inputsource = DSK6713_AIC23_INPUT_LINEIN;                // 0x011

#include <math.h>
#define PI 3.14159265358979

#define PTS 256                                                 // number of FFT points
#define RADIX 2                                                 // Division radix, FFT related
#define DELTA (2*PI)/PTS                                        // FFT frequency resolution

typedef struct Complex_tag {float real,imag;} COMPLEX;          // structure to represent complex numbers

#pragma DATA_ALIGN(W,sizeof(COMPLEX))                           // For FFT assembly code to work propoerly, the involved parameters should
#pragma DATA_ALIGN(outWin,sizeof(COMPLEX))                      // be aligned to 8byte memory locations. This tells compiler to exactly
#pragma DATA_ALIGN(h,sizeof(COMPLEX))                           // do that
#pragma DATA_ALIGN(bass,sizeof(COMPLEX))
#pragma DATA_ALIGN(mid1,sizeof(COMPLEX))
#pragma DATA_ALIGN(mid2,sizeof(COMPLEX))
#pragma DATA_ALIGN(treble,sizeof(COMPLEX))
COMPLEX W[PTS/RADIX] ;                                          // twiddle array, used in FFT
COMPLEX outWin[PTS];                                            // output window
COMPLEX h[PTS];                                                 // freq domain filter

#include "GraphicEQcoeff.h"                                     // time-domain FIR coefficients
COMPLEX bass[PTS], mid1[PTS], mid2[PTS], treble[PTS];                       // individual equalizer filters, frequency domain

float bass_gain = 1.0;                                          // Gains to be applied to each filter. Adjustable by DIP switches later
float mid1_gain = 1.0;
float mid2_gain = 1.0;
float treble_gain = 1.0;

//The two commented lines are needed of interrupt is used
//short buffercount = 0;                                            // buffer count for iobuffer
//short flag = 0;                                                   //set to indicate iobuffer full

const int s = 195;
float iobuffer[195];                                              // primary input/output buffer
float overlap[61];                                               // Vector to hold the tail
short i;                                                            // general index variable
float a, b;                                                         // Temp variables for complex multiply
short NUMCOEFFS = sizeof(lpcoeff)/sizeof(float);                    // number of filter coefficients in included filters
short iTwid[PTS/2]={0} ;                                            // Twiddle index array. Used in FFT



//  Don' forget to add following files to your project. Add only one of them!
//  Vectors_intr.asm -> for interrupt
//  Vectors_poll.asm -> for polling

//  You should implement your interrupt service routine (ISR) in here if your program uses interrupt.
/*
interrupt void c_int11(void)
{
    output_left_sample((short)(iobuffer[buffercount]));
    iobuffer[buffercount++] = (float)(input_left_sample());
    if (buffercount >= PTS/2) {
        buffercount = 0;
        flag = 1;
    }
}
*/


#include "digitrev_index.h"                                         // definition of digitrev_index function, used for FFT
#include "LED_fun.h"                                                // function to check dip switch status, turn LEDs on and off,
                                                                    // and set filter gains
//void main()
//{
//    for(i=0;i<PTS/2;i++){iTwid[i]=0;}
//    for(i=0;i<256-195;++i){overlap[i]=0;}
//
//    digitrev_index(iTwid, PTS/RADIX, RADIX);
//
//    for( i = 0; i < PTS/RADIX; i++ ) {
//        W[i].real = cos(DELTA*i);
//        W[i].imag = sin(DELTA*i);
//        }
//
//    bitrev(W, iTwid, PTS/RADIX);
//
//    for (i=0 ; i<PTS ; i++) {
//        bass[i].real = 0.0;     bass[i].imag = 0.0;
//        mid1[i].real = 0.0;     mid1[i].imag = 0.0;
//        mid2[i].real = 0.0;     mid2[i].imag = 0.0;
//        treble[i].real = 0.0;   treble[i].imag = 0.0;
//        }
//
//    for (i=0; i<NUMCOEFFS; i++)
//    {
//        bass[i].real = lpcoeff[i];
//        mid1[i].real = bp1coeff[i];
//        mid2[i].real = bp2coeff[i];
//        treble[i].real = hpcoeff[i];
//        }
//
//    cfftr2_dit(bass,W,PTS);
//    cfftr2_dit(mid1,W,PTS);
//    cfftr2_dit(mid2,W,PTS);
//    cfftr2_dit(treble,W,PTS);
//
//
//    DSK6713_LED_init();
//    DSK6713_DIP_init();
//
//    //comm_intr();
//    comm_poll();
//
//
//
//
//    while(1) {
//
//        LED_fun();
//
//        for(i=0;i<195;i++){
//            output_left_sample((short)(iobuffer[i]));
//            iobuffer[i] = (float)(input_left_sample());
//        }
//
//        // follwoing two lines are used for interrupt. Remember to comment the for loop above
//        //while (flag == 0);    //wait for iobuffer to fill up
//        //      flag = 0;
//
//
//        for (i=0 ; i<195 ; i++)   outWin[i].real = iobuffer[i];
//        for (i=195 ; i<PTS ; i++) outWin[i].real = 0;
//        for (i=0 ; i<PTS ; i++)     outWin[i].imag = 0.0;
//
//        cfftr2_dit(outWin,W,PTS);
//
//        for (i=0 ; i<PTS ; i++) {
//            h[i].real =     bass[i].real*bass_gain  +  mid1[i].real*mid1_gain
//                            +  mid2[i].real*mid2_gain + treble[i].real*treble_gain;
//            h[i].imag =     bass[i].imag*bass_gain  +  mid1[i].imag*mid1_gain
//                            +  mid2[i].imag*mid2_gain + treble[i].imag*treble_gain;
//        }
//
//        for (i=0; i<PTS; i++) {
//            a = outWin[i].real;
//            b = outWin[i].imag;
//            outWin[i].real = h[i].real*a - h[i].imag*b;
//            outWin[i].imag = h[i].real*b + h[i].imag*a;
//        }
//
//        icfftr2_dif(outWin,W,PTS);
//
//        for (i=0 ; i<PTS ; i++) outWin[i].real /= PTS;
//
//        // Overlap Add Method
//
//        for(i=0;i<61;i++){
//            iobuffer[i]=outWin[i].real+overlap[i];
//            overlap[i]=outWin[i+61].real;
//        }
//        for(i=61;i<195;++i)
//            iobuffer[i] = outWin[i].real;
//    }
//
//}



// Overlap Save
void main()
{
    for(i=0;i<PTS/2;i++){iTwid[i]=0;}
    for(i=0;i<256-195;++i){overlap[i]=0;}

    digitrev_index(iTwid, PTS/RADIX, RADIX);

    for( i = 0; i < PTS/RADIX; i++ ) {
        W[i].real = cos(DELTA*i);
        W[i].imag = sin(DELTA*i);
        }

    bitrev(W, iTwid, PTS/RADIX);

    for (i=0 ; i<PTS ; i++) {
        bass[i].real = 0.0;     bass[i].imag = 0.0;
        mid1[i].real = 0.0;     mid1[i].imag = 0.0;
        mid2[i].real = 0.0;     mid2[i].imag = 0.0;
        treble[i].real = 0.0;   treble[i].imag = 0.0;
        }

    for (i=0; i<NUMCOEFFS; i++)
    {
        bass[i].real = lpcoeff[i];
        mid1[i].real = bp1coeff[i];
        mid2[i].real = bp2coeff[i];
        treble[i].real = hpcoeff[i];
        }

    cfftr2_dit(bass,W,PTS);
    cfftr2_dit(mid1,W,PTS);
    cfftr2_dit(mid2,W,PTS);
    cfftr2_dit(treble,W,PTS);


    DSK6713_LED_init();
    DSK6713_DIP_init();

    //comm_intr();
    comm_poll();




    while(1) {

        LED_fun();

        for(i=0;i<195;i++){
            output_left_sample((short)(iobuffer[i]));
            iobuffer[i] = (float)(input_left_sample());
        }

        // follwoing two lines are used for interrupt. Remember to comment the for loop above
        //while (flag == 0);    //wait for iobuffer to fill up
        //      flag = 0;

        for(i = 0; i< 61; i++) outWin[i].real = overlap[i];
        for (i=61 ; i<PTS ; i++)   outWin[i].real = iobuffer[i-61];
//        for (i=195 ; i<PTS ; i++) outWin[i].real = 0;
        for (i=0 ; i<PTS ; i++)     outWin[i].imag = 0.0;

        cfftr2_dit(outWin,W,PTS);

        for (i=0 ; i<PTS ; i++) {
            h[i].real =     bass[i].real*bass_gain  +  mid1[i].real*mid1_gain
                            +  mid2[i].real*mid2_gain + treble[i].real*treble_gain;
            h[i].imag =     bass[i].imag*bass_gain  +  mid1[i].imag*mid1_gain
                            +  mid2[i].imag*mid2_gain + treble[i].imag*treble_gain;
        }

        for (i=0; i<195; i++) {
            a = outWin[i].real;
            b = outWin[i].imag;
            outWin[i].real = h[i].real*a - h[i].imag*b;
            outWin[i].imag = h[i].real*b + h[i].imag*a;
        }

        icfftr2_dif(outWin,W,PTS);

        for (i=0 ; i<PTS ; i++) outWin[i].real /= PTS;


        for (i=0; i < 61; ++i)
        	overlap[i] = iobuffer[195-61+i];
        for(i = 0; i < 195; ++i)
        	iobuffer[i] = outWin[i+61].real;


    }

}
