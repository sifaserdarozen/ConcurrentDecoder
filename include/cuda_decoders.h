#ifndef _CUDA_DECODER_LIB_H
#define _CUDA_DECODER_LIB_H

#define THREAD_PER_BLOCK 256

void CudaGpuInitialize();
int CudaG711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int CudaG711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);
int CudaG722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr);

#endif