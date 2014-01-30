#include <iostream>
#include <cuda.h>
#include "cuda_decoders.h"

#include "g722decoder.h"

#define checkCudaErrors(val) CheckErrors( (val), __FILE__, __LINE__)

void CheckErrors(cudaError_t cuda_error, const char* const file, const int line)
{
    if (cuda_error != cudaSuccess)
	{
        std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
        std::cerr << cudaGetErrorString(cuda_error) << std::endl;
        exit(1);
    }		
}

__global__ void CudaKernelG711aToPcm(unsigned char* d_alaw_data_ptr, short int* d_pcm_data_ptr)
{
    unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);
	
	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;

	for (int k=0; k<160; k++)
	{
	    alaw_data = d_alaw_data_ptr[idx+k];	
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
		
		d_pcm_data_ptr[idx+k]=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

// memory coalesced version of alaw to pcm conversion
__global__ void CudaKernelG711aToPcmCM(unsigned char* d_alaw_data_ptr, short int* d_pcm_data_ptr)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;
	
	short int quantization_value;
	short int quantization_segment;
	unsigned char alaw_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
	    alaw_data = d_alaw_data_ptr[idx];	
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
		
		d_pcm_data_ptr[idx]=((alaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

__global__ void CudaKernelG711uToPcm(unsigned char* d_ulaw_data_ptr, short int* d_pcm_data_ptr)
{
    unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);
	
	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;

	for (int k=0; k<160; k++)
	{
		ulaw_data=~(d_ulaw_data_ptr[idx+k]);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);
	
		d_pcm_data_ptr[idx+k]=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

__global__ void CudaKernelG711uToPcmCM(unsigned char* d_ulaw_data_ptr, short int* d_pcm_data_ptr)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;
		
	short int quantization_value;
	short int quantization_segment;
	unsigned char ulaw_data;

	for (int k=0; k<160; k++, idx+=total_threads)
	{
		ulaw_data=~(d_ulaw_data_ptr[idx]);

		quantization_value= (ulaw_data & (0xf)) << 4;
		quantization_segment = ((unsigned)ulaw_data & (0x70)) >> (4);

		quantization_value += 0x0084;
		quantization_value <<= quantization_segment;

		quantization_value-=(32);
	
		d_pcm_data_ptr[idx]=((ulaw_data & (0x80))?quantization_value:-quantization_value);
	}
}

__global__ void CudaKernelG722ToPcm(unsigned char* d_g722_data_ptr, short int* d_pcm_data_ptr)
{
    unsigned int idx = 160*(threadIdx.x + blockDim.x * blockIdx.x);
	
	unsigned char g722_data;

	for (int k=0; k<160; k++)
	{
		g722_data=d_g722_data_ptr[idx+k];

		d_pcm_data_ptr[idx+k]=0;
	}
}

__device__ short int CudaConvertLongToShort(int in_value)
{
	if (in_value > 32767)
		return 32767;
	else if (in_value < -32768)
		return -32768;
	else
		return (short)in_value;
}

__global__ void CudaKernelG722ToPcmCM(unsigned char* d_g722_data_ptr, short int* d_pcm_data_ptr, int* d_band_data_ptr, int* d_g722_consts_ptr, unsigned int no_of_data)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int total_threads = blockDim.x * gridDim.x;
		
	unsigned char g722_data;
	
	int number_of_chunks = no_of_data/160;
	if (idx >= number_of_chunks)
	    return;
		
	// pointers for constants, maybe copy these to shared mem
	int* wl = d_g722_consts_ptr;
	int* rl42 = d_g722_consts_ptr + 8;
	int* ilb = rl42 + 16;
	int* qm4 = ilb + 32;
	int* qm6 = qm4 + 16;

	
	// copy band data to local variables
	int band_s = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_sp = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_sz = d_band_data_ptr[idx];
	idx += number_of_chunks;
	
	int band_r[3], band_a[3], band_ap[3], band_p[3];
	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
	    band_r[k] = d_band_data_ptr[idx];
		band_a[k] = d_band_data_ptr[idx + 3*number_of_chunks];
		band_ap[k] = d_band_data_ptr[idx + 6*number_of_chunks];
	    band_p[k] = d_band_data_ptr[idx + 9*number_of_chunks];
	}
	idx += 9*number_of_chunks;
	
	int band_d[7], band_b[7], band_bp[7], band_sg[7];
	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
	    band_d[k] = d_band_data_ptr[idx];
		band_b[k] = d_band_data_ptr[idx + 7*number_of_chunks];
		band_bp[k] = d_band_data_ptr[idx + 14*number_of_chunks];
	    band_sg[k] = d_band_data_ptr[idx + 21*number_of_chunks];
	}
	idx += 21*number_of_chunks;
	
	int band_nb = d_band_data_ptr[idx];
	idx += number_of_chunks;
	int band_det = d_band_data_ptr[idx];
	//band_det=32;
	
	int dlowt;
	int rlow;
	int wd1;
	int wd2;
	int wd3;
	
	idx = threadIdx.x + blockDim.x * blockIdx.x;
	for (int k=0; k<160; k++, idx+=number_of_chunks)
	{
		g722_data=d_g722_data_ptr[idx];

	    wd1 = g722_data & 0x3F;
	    wd2 = qm6[wd1];
	    wd1 >>= 2;

	    /********************** Block 5 *******************/
	    // INVQBL (ITU page 43), compute quantized difference signal for the decoder output in the lower sub-band
	    wd2 = (band_det * wd2) >> 15;
	    // RECONS ( ITU page 41), compute reconstructed signal for the adaptive predictor
	    rlow = band_s + wd2;
	
	    /********************** Block 6 ********************/
	    // LIMIT (ITU page 44), limit the output reconstructed signal
	    if (rlow > 16383)
		    rlow = 16383;
	    else if (rlow < -16384)
		    rlow = -16384;

	    /********************** Block 2 ***********************/	
	    // INVQAL (ITU page 37), compute the quantized differences signal for the adaptive predictor in the lower sub-band
	    wd2 = qm4[wd1];
	    dlowt = (band_det * wd2) >> 15;

	    /********************** Block 3 ************************/
	    // LOGSCL (ITU page 38), update the logarithmic quantizer scale factor in the lower sub-band
	    wd2 = rl42[wd1];
	    wd1 = (band_nb * 127) >> 7;
	    wd1 += wl[wd2];
	    if (wd1 < 0)
		    wd1 = 0;
	    else if (wd1 > 18432)
		    wd1 = 18432;
	    band_nb = wd1;

	    // SCALEL (ITU page 38), compute the quantizer scale factor in the lower sub-band 
	    wd1 = (band_nb >> 6) & 31;
	    wd2 = 8 - (band_nb >> 11);
	    wd3 = (wd2 < 0)	 ?  (ilb[wd1] << -wd2)	:  (ilb[wd1] >> wd2);
	    band_det = wd3 << 2;

	    /********************** Block 4 **************************/

	    // RECONS (ITU page 41), compute reconstructed signal for the adaptive predictor
	    band_d[0] = dlowt;
	    band_r[0] = CudaConvertLongToShort(band_s + dlowt);

	    // PARREC (ITU page 40), compute partially reconstructed signal
	    band_p[0] = CudaConvertLongToShort(band_sz + dlowt);

	    // UPPOL2 (ITU page 41), update second predictor coefficient
	    int i;  // loop variable
	    for (i = 0;	 i < 3;	 i++)
		    band_sg[i] = band_p[i] >> 15;
	    wd1 = CudaConvertLongToShort(band_a[1] << 2);

	    wd2 = (band_sg[0] == band_sg[1])	?  -wd1	 :  wd1;
	    if (wd2 > 32767)
		    wd2 = 32767;
	    wd3 = (band_sg[0] == band_sg[2])	?  128	:  -128;
	    wd3 += (wd2 >> 7);
	    wd3 += (band_a[2]*32512) >> 15;
	    if (wd3 > 12288)
	    	wd3 = 12288;
	    else if (wd3 < -12288)
		    wd3 = -12288;
	    band_ap[2] = wd3;

	    // UPPOL1 (ITU page 42), update first predictor coefficient
	    band_sg[0] = band_p[0] >> 15;
	    band_sg[1] = band_p[1] >> 15;
	    wd1 = (band_sg[0] == band_sg[1])	?  192	:  -192;
	    wd2 = (band_a[1]*32640) >> 15;

	    band_ap[1] = CudaConvertLongToShort(wd1 + wd2);
	    wd3 = CudaConvertLongToShort(15360 - band_ap[2]);
	    if (band_ap[1] > wd3)
		    band_ap[1] = wd3;
	    else if (band_ap[1] < -wd3)
		    band_ap[1] = -wd3;

	    // UPZERO (ITU page 41), update sixth order predictor coefficients
	    wd1 = (dlowt == 0)  ?  0  :  128;
	    band_sg[0] = dlowt >> 15;
	    for (i = 1;	 i < 7;	 i++)
	    {
		    band_sg[i] = band_d[i] >> 15;
		    wd2 = (band_sg[i] == band_sg[0])  ?  wd1  :  -wd1;
		    wd3 = (band_b[i]*32640) >> 15;
		    band_bp[i] = CudaConvertLongToShort(wd2 + wd3);
	    }

	    // DELAYA (ITU page 38), memory block delay 
	    for (i = 6;	 i > 0;	 i--)
	    {
		    band_d[i] = band_d[i - 1];
		    band_b[i] = band_bp[i];
	    }

	    for (i = 2;	 i > 0;	 i--)
	    {
		    band_r[i] = band_r[i - 1];
		    band_p[i] = band_p[i - 1];
		    band_a[i] = band_ap[i];
	    }

	    // FILTEP (ITU page 43), compute predictor output signal, poles
	    wd1 = CudaConvertLongToShort(band_r[1] + band_r[1]);
	    wd1 = (band_a[1]*wd1) >> 15;
	    wd2 = CudaConvertLongToShort(band_r[2] + band_r[2]);
	    wd2 = (band_a[2]*wd2) >> 15;
	    band_sp = CudaConvertLongToShort(wd1 + wd2);

	    // FILTEZ (ITU page 42), compute predictor output signal, zeros
	    band_sz = 0;
	    for (i = 6;	 i > 0;	 i--)
	    {
		    wd1 = CudaConvertLongToShort(band_d[i] + band_d[i]);
		    band_sz += (band_b[i]*wd1) >> 15;
	    }
	    band_sz = CudaConvertLongToShort(band_sz);

	    // PREDIC (ITU page 43), compute predictor output value
	    band_s = CudaConvertLongToShort(band_sp + band_sz);

		d_pcm_data_ptr[idx]=(short int)rlow;
	}
	
	
	
	// copy local variables back to band data
	idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	d_band_data_ptr[idx] = band_s;
	idx += number_of_chunks;
	band_sp = d_band_data_ptr[idx]= band_sp;
	idx += number_of_chunks;
	d_band_data_ptr[idx] = band_sz;
	idx += number_of_chunks;
	
	for(int k=0; k<3; k++, idx+=number_of_chunks)
	{
	    d_band_data_ptr[idx] = band_r[k];
		d_band_data_ptr[idx + 3*number_of_chunks] = band_a[k];
		d_band_data_ptr[idx + 6*number_of_chunks] = band_ap[k];
	    d_band_data_ptr[idx + 9*number_of_chunks] = band_p[k];
	}
	idx += 9*number_of_chunks;
	
	for(int k=0; k<7; k++, idx+=number_of_chunks)
	{
	    d_band_data_ptr[idx] = band_d[k];
		d_band_data_ptr[idx + 7*number_of_chunks] = band_b[k];
		d_band_data_ptr[idx + 14*number_of_chunks] = band_bp[k];
	    d_band_data_ptr[idx + 21*number_of_chunks] = band_sg[k];
	}
	idx += 21*number_of_chunks;
	
	d_band_data_ptr[idx] = band_nb;
	idx += number_of_chunks;
	d_band_data_ptr[idx] = band_det;
	idx += number_of_chunks;
}

int CudaG711aToPcm(const unsigned char* alaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
    dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);
	
	unsigned int size_of_alaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);
	
	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_alaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);
	
	unsigned char* d_alaw_data_ptr = NULL;
    cudaMalloc((void**)&d_alaw_data_ptr, size_of_d_alaw_data);
	checkCudaErrors(cudaGetLastError());
	
	short int* d_pcm_data_ptr = NULL;
    cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());
	
    cudaMemcpy(d_alaw_data_ptr, alaw_data_ptr, size_of_alaw_data, cudaMemcpyHostToDevice);
    checkCudaErrors(cudaGetLastError());
	
    // launch kernel here
	CudaKernelG711aToPcmCM <<< grid_dim, block_dim >>> (d_alaw_data_ptr, d_pcm_data_ptr);
    checkCudaErrors(cudaGetLastError());
	
    cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
    checkCudaErrors(cudaGetLastError());
		
    cudaFree(d_alaw_data_ptr);
	checkCudaErrors(cudaGetLastError());
    
	cudaFree(d_pcm_data_ptr);
    checkCudaErrors(cudaGetLastError());	

    return 0;
}

int CudaG711uToPcm(const unsigned char* ulaw_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
    dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);
	
	unsigned int size_of_ulaw_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);
	
	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_ulaw_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);
	
	unsigned char* d_ulaw_data_ptr = NULL;
    cudaMalloc((void**)&d_ulaw_data_ptr, size_of_d_ulaw_data);
	checkCudaErrors(cudaGetLastError());
	
	short int* d_pcm_data_ptr = NULL;
    cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());
	
    cudaMemcpy(d_ulaw_data_ptr, ulaw_data_ptr, size_of_ulaw_data, cudaMemcpyHostToDevice);
    checkCudaErrors(cudaGetLastError());
	
    // launch kernel here
	CudaKernelG711uToPcmCM <<< grid_dim, block_dim >>> (d_ulaw_data_ptr, d_pcm_data_ptr);
    checkCudaErrors(cudaGetLastError());
	
    cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
    checkCudaErrors(cudaGetLastError());
		
    cudaFree(d_ulaw_data_ptr);
	checkCudaErrors(cudaGetLastError());
    
	cudaFree(d_pcm_data_ptr);
    checkCudaErrors(cudaGetLastError());	

    return 0;
}

int CudaG722ToPcm(const unsigned char* g722_data_ptr, int* band_data_ptr, unsigned int no_of_data, short int* pcm_data_ptr)
{
    dim3 block_dim(THREAD_PER_BLOCK, 1, 1);
	dim3 grid_dim(ceil(no_of_data/((float)(block_dim.x * 160))), 1, 1);
	
	unsigned int size_of_g722_data = no_of_data * sizeof(unsigned char);
	unsigned int size_of_pcm_data = no_of_data * sizeof(short int);
	
	unsigned int no_of_d_data = grid_dim.x * block_dim.x * 160;
	unsigned int size_of_d_g722_data = no_of_d_data * sizeof(unsigned char);
	unsigned int size_of_d_pcm_data = no_of_d_data * sizeof(short int);
	
	unsigned char* d_g722_data_ptr = NULL;
    cudaMalloc((void**)&d_g722_data_ptr, size_of_d_g722_data);
	checkCudaErrors(cudaGetLastError());
	
	short int* d_pcm_data_ptr = NULL;
    cudaMalloc((void**)&d_pcm_data_ptr, size_of_d_pcm_data);
	checkCudaErrors(cudaGetLastError());
	
    cudaMemcpy(d_g722_data_ptr, g722_data_ptr, size_of_g722_data, cudaMemcpyHostToDevice);
    checkCudaErrors(cudaGetLastError());
	
	// calculate space for band, 45 integers per thread
	unsigned int number_of_d_band_data = grid_dim.x * block_dim.x;
	unsigned int size_of_d_band_data = number_of_d_band_data * 45 * sizeof(int);
	unsigned int number_of_band_data = no_of_data/160;
	unsigned int size_of_band_data = number_of_band_data * 45 * sizeof(int);
	
	//std::cout << "size of band data " << size_of_d_band_data << std::endl;
	
	int* d_band_data_ptr = NULL;
    cudaMalloc((void**)&d_band_data_ptr, size_of_d_band_data);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_band_data_ptr, band_data_ptr, size_of_band_data, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	
	
	unsigned int size_of_d_g722_consts = sizeof(g722_consts);
	int* d_g722_consts_ptr = NULL;
    cudaMalloc((void**)&d_g722_consts_ptr, size_of_d_g722_consts);
	checkCudaErrors(cudaGetLastError());
	cudaMemcpy(d_g722_consts_ptr, g722_consts, size_of_d_g722_consts, cudaMemcpyHostToDevice);
	checkCudaErrors(cudaGetLastError());	
 
	
    // launch kernel here
	CudaKernelG722ToPcmCM <<< grid_dim, block_dim >>> (d_g722_data_ptr, d_pcm_data_ptr, d_band_data_ptr, d_g722_consts_ptr, no_of_data);
    checkCudaErrors(cudaGetLastError());
	
	//std::cout << "host pcm data size   : " << size_of_pcm_data << std::endl;
	//std::cout << "device pcm data size : " << size_of_d_pcm_data << std::endl;
	//std::cout << "no of data           : " << no_of_data << std::endl;
	//std::cout << "number of threads    : " << no_of_d_data << std::endl;
	
    cudaMemcpy(pcm_data_ptr, d_pcm_data_ptr, size_of_pcm_data, cudaMemcpyDeviceToHost);
    checkCudaErrors(cudaGetLastError());
	
	cudaMemcpy(band_data_ptr, d_band_data_ptr, size_of_band_data, cudaMemcpyDeviceToHost);
	checkCudaErrors(cudaGetLastError());
	
	cudaFree(d_g722_consts_ptr);
	checkCudaErrors(cudaGetLastError());
	
	cudaFree(d_band_data_ptr);
	checkCudaErrors(cudaGetLastError());
		
    cudaFree(d_g722_data_ptr);
	checkCudaErrors(cudaGetLastError());
    
	cudaFree(d_pcm_data_ptr);
    checkCudaErrors(cudaGetLastError());	

    return 0;
}

void CudaGpuInitialize()
{
    unsigned int size_of_d_dummy_data = 1000000;
    unsigned char* d_dummy_data_ptr = NULL;

    //unsigned int device_count = cudaGetDeviceCount();
    //checkCudaErrors(cudaGetLastError()); 
    
    //std::cout << "device count : " << device_count << std::endl
	
    cudaMalloc((void**)&d_dummy_data_ptr, size_of_d_dummy_data);
    checkCudaErrors(cudaGetLastError());
	
    cudaMemset((void*)d_dummy_data_ptr, 0, size_of_d_dummy_data);
    checkCudaErrors(cudaGetLastError());
	
    cudaFree(d_dummy_data_ptr);
    checkCudaErrors(cudaGetLastError());
}
