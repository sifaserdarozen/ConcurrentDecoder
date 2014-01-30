#include "cpu_decoders.h"
#include "g722decoder.h"
#include <iostream>

int G711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!alaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	//double normalizing_ratio=((double)(0x9ffff))/((double)(0x1ffff)); // this will map 13 bit pcm to pseudo 16 bit pcm
	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;
	short int pcm_data;
	for (int k=0; k<no_of_data; k++)
	{
		alaw_data=*alaw_data_ptr++;
		alaw_data^=0x55;

		quantization_value= (alaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)alaw_data & (0x70)) >> (4);
		switch (quantization_segment)
		{
		case 0: 
			quantization_value+=(0x0008);
			break;
		case 1:
			quantization_value+=(0x0108);
			break;
		default:
			quantization_value+=(0x0108);
			quantization_value <<= (quantization_segment-1);
		};
		
		//*pcm_vector++=((alaw_data & (0x80))?quantization_value:-quantization_value)*normalizing_ratio;
		*pcm_data_ptr++=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}

	return 0;
}

int PcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!alaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << " at int PcmToAlaw()" << " at : int CallType::PcmToAlaw" << std::endl;
		return -1;
	}

	short int quantization_segment;
	short int quantization_value;

	short int pcm_data;
	unsigned char alaw_data;

	for(int k=0; k<no_of_data; k++)
	{
		pcm_data=*pcm_data_ptr++;
		quantization_value=(pcm_data<0) ? ((~pcm_data)>>4) : (pcm_data>>4);
		
		if(quantization_value>15)
		{
			quantization_segment=1;
			while(quantization_value>(16+15))
			{
				quantization_value>>=1;
				quantization_segment++;
			}
			quantization_value-=16;

			alaw_data=quantization_value + (quantization_segment << 4);
		}

		if(pcm_data>=0)
			alaw_data |= 0x80;

		alaw_data^=0x55;

		*alaw_data_ptr++=alaw_data;
	}

	return 0;
}

int G711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!ulaw_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}

	//double normalizing_ratio=((double)(0x9ffff))/((double)(0x1ffff)); // this will map 13 bit pcm to pseudo 16 bit pcm
	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;
	short int pcm_data;
	for (int k=0; k<no_of_data; k++)
	{
		ulaw_data=~(*ulaw_data_ptr++);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);
	
		//*pcm_vector++=((alaw_data & (0x80))?quantization_value:-quantization_value)*normalizing_ratio;
		*pcm_data_ptr++=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}

	return 0;
}

int G722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
	if (no_of_data < 0)
	{
		std::cerr << "Number of elements is : " << no_of_data << std::endl;
		return -1;
	}

	if (!g722_data_ptr)
	{
		std::cerr << "alaw vector is null" << std::endl;
		return -1;
	}

	if (!pcm_data_ptr)
	{
		std::cerr << "pcm vector is null" << std::endl;
		return -1;
	}
	
	int number_of_chunks = no_of_data/160;
	
	for(unsigned int k=0; k<number_of_chunks; k++)
	{
	
	    unsigned char g722_data_arr[160];
		unsigned short pcm_data_arr[160];
		int band_data[45];
		
		for(int j=0; j<160; j++)
		    g722_data_arr[j] = g722_data_ptr[k + number_of_chunks*j];
		
        for(int j=0; j<45; j++)
             band_data[j]= band_data_ptr[k + number_of_chunks*j];		
			
		G722DecoderType g722_decoder;
		g722_decoder.Decode(g722_data_arr, band_data, pcm_data_arr, 160);
	
	    for(int j=0; j<160; j++)
		    pcm_data_ptr[k + number_of_chunks*j] = pcm_data_arr[j];

        for(int j=0; j<45; j++)
            band_data_ptr[k+ number_of_chunks*j] = band_data[j];		
	}
	
    return 0;
}
