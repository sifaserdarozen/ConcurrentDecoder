#include "g722decoder.h"
#include <iostream>
#include <cstring>    // needed for g++ for std::memcpy

G722DecoderType::G722DecoderType()
{
    Reset();
}

G722DecoderType::~G722DecoderType()
{

}

short int G722DecoderType::ConvertLongToShort(long int in_value)
{
	if (in_value > 32767)
		return 32767;
	else if (in_value < -32768)
		return -32768;
	else
		return (short int)in_value;
}

void G722DecoderType::Decode(unsigned char* g722_data_ptr, int* band_data, unsigned short* pcm_data_ptr, unsigned int no_of_data)
{
    std::memcpy((void*)&band, (void*)band_data, 45*sizeof(int));

	int dlowt;
	int rlow;
//	int ihigh;
	int wd1;
	int wd2;
	int wd3;
	unsigned char g722_data;
	
	for (unsigned int index=0; index<no_of_data; index++)
	{
	    g722_data = g722_data_ptr[index];

	    wd1 = g722_data & 0x3F;
//	    ihigh = (g722_data >> 6) & 0x03;
	    wd2 = qm6[wd1];
	    wd1 >>= 2;

	    /********************** Block 5 *******************/
	    // INVQBL (ITU page 43), compute quantized difference signal for the decoder output in the lower sub-band
	    wd2 = (band.det * wd2) >> 15;
	    // RECONS ( ITU page 41), compute reconstructed signal for the adaptive predictor
	    rlow = band.s + wd2;
	
	    /********************** Block 6 ********************/
	    // LIMIT (ITU page 44), limit the output reconstructed signal
	    if (rlow > 16383)
		    rlow = 16383;
	    else if (rlow < -16384)
		    rlow = -16384;

	    /********************** Block 2 ***********************/	
	    // INVQAL (ITU page 37), compute the quantized differences signal for the adaptive predictor in the lower sub-band
	    wd2 = qm4[wd1];
	    dlowt = (band.det * wd2) >> 15;

	    /********************** Block 3 ************************/
	    // LOGSCL (ITU page 38), update the logarithmic quantizer scale factor in the lower sub-band
	    wd2 = rl42[wd1];
	    wd1 = (band.nb * 127) >> 7;
	    wd1 += wl[wd2];
	    if (wd1 < 0)
		    wd1 = 0;
	    else if (wd1 > 18432)
		    wd1 = 18432;
	    band.nb = wd1;

	    // SCALEL (ITU page 38), compute the quantizer scale factor in the lower sub-band 
	    wd1 = (band.nb >> 6) & 31;
	    wd2 = 8 - (band.nb >> 11);
	    wd3 = (wd2 < 0)	 ?  (ilb[wd1] << -wd2)	:  (ilb[wd1] >> wd2);
	    band.det = wd3 << 2;

	    /********************** Block 4 **************************/

	    // RECONS (ITU page 41), compute reconstructed signal for the adaptive predictor
	    band.d[0] = dlowt;
	    band.r[0] = ConvertLongToShort(band.s + dlowt);

	    // PARREC (ITU page 40), compute partially reconstructed signal
	    band.p[0] = ConvertLongToShort(band.sz + dlowt);

	    // UPPOL2 (ITU page 41), update second predictor coefficient
	    int i;  // loop variable
	    for (i = 0;	 i < 3;	 i++)
		    band.sg[i] = band.p[i] >> 15;
	    wd1 = ConvertLongToShort(band.a[1] << 2);

	    wd2 = (band.sg[0] == band.sg[1])	?  -wd1	 :  wd1;
	    if (wd2 > 32767)
		    wd2 = 32767;
	    wd3 = (band.sg[0] == band.sg[2])	?  128	:  -128;
	    wd3 += (wd2 >> 7);
	    wd3 += (band.a[2]*32512) >> 15;
	    if (wd3 > 12288)
	    	wd3 = 12288;
	    else if (wd3 < -12288)
		    wd3 = -12288;
	    band.ap[2] = wd3;

	    // UPPOL1 (ITU page 42), update first predictor coefficient
	    band.sg[0] = band.p[0] >> 15;
	    band.sg[1] = band.p[1] >> 15;
	    wd1 = (band.sg[0] == band.sg[1])	?  192	:  -192;
	    wd2 = (band.a[1]*32640) >> 15;

	    band.ap[1] = ConvertLongToShort(wd1 + wd2);
	    wd3 = ConvertLongToShort(15360 - band.ap[2]);
	    if (band.ap[1] > wd3)
		    band.ap[1] = wd3;
	    else if (band.ap[1] < -wd3)
		    band.ap[1] = -wd3;

	    // UPZERO (ITU page 41), update sixth order predictor coefficients
	    wd1 = (dlowt == 0)  ?  0  :  128;
	    band.sg[0] = dlowt >> 15;
	    for (i = 1;	 i < 7;	 i++)
	    {
		    band.sg[i] = band.d[i] >> 15;
		    wd2 = (band.sg[i] == band.sg[0])  ?  wd1  :  -wd1;
		    wd3 = (band.b[i]*32640) >> 15;
		    band.bp[i] = ConvertLongToShort(wd2 + wd3);
	    }

	    // DELAYA (ITU page 38), memory block delay 
	    for (i = 6;	 i > 0;	 i--)
	    {
		    band.d[i] = band.d[i - 1];
		    band.b[i] = band.bp[i];
	    }

	    for (i = 2;	 i > 0;	 i--)
	    {
		    band.r[i] = band.r[i - 1];
		    band.p[i] = band.p[i - 1];
		    band.a[i] = band.ap[i];
	    }

	    // FILTEP (ITU page 43), compute predictor output signal, poles
	    wd1 = ConvertLongToShort(band.r[1] + band.r[1]);
	    wd1 = (band.a[1]*wd1) >> 15;
	    wd2 = ConvertLongToShort(band.r[2] + band.r[2]);
	    wd2 = (band.a[2]*wd2) >> 15;
	    band.sp = ConvertLongToShort(wd1 + wd2);

	    // FILTEZ (ITU page 42), compute predictor output signal, zeros
	    band.sz = 0;
	    for (i = 6;	 i > 0;	 i--)
	    {
		    wd1 = ConvertLongToShort(band.d[i] + band.d[i]);
		    band.sz += (band.b[i]*wd1) >> 15;
	    }
	    band.sz = ConvertLongToShort(band.sz);

	    // PREDIC (ITU page 43), compute predictor output value
	    band.s = ConvertLongToShort(band.sp + band.sz);

	    pcm_data_ptr[index]= (short int)rlow;
	}
	std::memcpy((void*)band_data, (void*)&band, 45*sizeof(int));
}

void G722DecoderType::Reset()
{
	band.s = 0;
	band.sp = 0;
	band.sz = 0;
	
	for (int k=0; k<3; k++)
	{
	    band.r[k] = 0;
	    band.a[k] = 0;
	    band.ap[k] = 0;
	    band.p[k] = 0;
	}
	
	for (int k=0; k<7; k++)
	{
	    band.d[k] = 0;
		band.b[k] = 0;
		band.bp[k] = 0;
		band.sg[k] = 0;
	}
	
	band.nb = 0;
    band.det = 32;
}
