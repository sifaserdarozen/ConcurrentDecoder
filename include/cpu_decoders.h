#ifndef _CPU_DECODER_LIB_H
#define _CPU_DECODER_LIB_H

int G711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int G711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int G722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);

int PcmToG711a(const short int* pcm_data_ptr, unsigned int no_of_data, unsigned char* alaw_data_ptr);

#endif