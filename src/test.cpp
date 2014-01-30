#include <iostream>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <cstring>    // for g++ std::memcpy
#include <cmath>
#include <iomanip>    // for std::setprecision
#include <vector>
#include <fstream>

#include "sample_data.h"
#include "cuda_decoders.h"
#include "cpu_decoders.h"

#define NO_OF_TEST_PACKET 1000000

#define DEFAULT_START_OF_ITERATIONS 250
#define DEFAULT_END_OF_ITERATIONS 128000
#define DEFAULT_NUMBER_OF_ITERATIONS 10
#define DEFAULT_NUMBER_OF_AVERAGING 1000

enum DataType {
    NAIVE,
    COALESCED
};

int GenerateG711aData(unsigned char* alaw_data_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!alaw_data_ptr)
        return -1;
		
    unsigned char* dst_ptr = alaw_data_ptr;
    const unsigned char* src_ptr = NULL;
    
    unsigned int no_of_source_data = sizeof(sample_g711a_data) / sizeof(sample_g711a_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711a_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
//          *dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;
        } 	
	
        // copy remeinder datas
        src_ptr = sample_g711a_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
//          *dst_ptr++ = *src_ptr++;	
            *dst_ptr++ = rand() % 256;
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g711u_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
    	
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }				
    }
			
    return 0;
}

int GenerateG711uData(unsigned char* ulaw_data_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!ulaw_data_ptr)
        return -1;
		
    unsigned char* dst_ptr = ulaw_data_ptr;
    const unsigned char* src_ptr = NULL;
    
    unsigned int no_of_source_data = sizeof(sample_g711u_data) / sizeof(sample_g711u_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711u_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
        }
	
        // copy remeinder datas
        src_ptr = sample_g711a_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
            //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;			
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g711u_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }				
    }
	
    return 0;
}

int GenerateG722Data(unsigned char* g722_data_ptr, int* g722_band_ptr, unsigned int no_of_data, DataType choice = NAIVE)
{
    if(!g722_data_ptr)
        return -1;

    unsigned char* dst_ptr = g722_data_ptr;
    const unsigned char* src_ptr = NULL;

    unsigned int no_of_source_data = sizeof(sample_g722_data) / sizeof(sample_g722_data[0]);
    unsigned int no_of_full_repetition = no_of_data/no_of_source_data;
    unsigned int last_repetition_remeinder = no_of_data%no_of_source_data;
    unsigned int no_of_total_repetitions = ceil(no_of_data/((float)no_of_source_data));
    unsigned int band_size_of_total_repetitions = 45 * sizeof(int) * no_of_total_repetitions;

    std::memset(g722_band_ptr, 0, band_size_of_total_repetitions);

    if(NAIVE == choice)
    {
        // copy full repetitions
        for(unsigned int k=0; k<no_of_full_repetition; k++)
        {
            src_ptr = sample_g711u_data;
            for(unsigned int j=0; j<no_of_source_data; j++)
                //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;
        }	

        // copy remeinder datas
        src_ptr = sample_g722_data;
        for(unsigned int j=0; j<last_repetition_remeinder; j++)
            //*dst_ptr++ = *src_ptr++;
            *dst_ptr++ = rand() % 256;

        for(unsigned int k=0; k<no_of_total_repetitions; k++)
            g722_band_ptr[k*45 + 44] = 32;		
    }
    else if (COALESCED == choice)
    {
        // copy full repetitions
        src_ptr = sample_g722_data;
        for(unsigned int j=0; j<no_of_source_data; j++)
        {
            for(unsigned int k=0; k<no_of_full_repetition; k++)
            {
                //*dst_ptr++ = *src_ptr;
                *dst_ptr++ = rand() % 256;				
            }
            if(j<last_repetition_remeinder)
            {
                //*dst_ptr++ = *src_ptr++;
                *dst_ptr++ = rand() % 256;
            }
            else
            {
                src_ptr++;
            }
        }

         for(unsigned int k=0; k<no_of_total_repetitions; k++)
            g722_band_ptr[k + 44*no_of_total_repetitions] = 32;				
    }
	
    return 0;
}

template <class T>
bool CompareArrays(const T* first_data_ptr, const T* second_data_ptr, unsigned int no_of_data)
{
    for(unsigned int k=0; k<no_of_data; k++)
    {
        if(*first_data_ptr++ != *second_data_ptr++)
        {
            std::cout << "match failed at index " << k << std::endl;
            return false;
        }			
    }
    
    std::cout << "match succedded..." << std::endl;
    return true;
}

int main(int argc, char *argv[])
{
    std::vector<int> iteration_vector;
    std::vector<int>::iterator itr;

    std::vector<float> cpu_alaw_vector;
    std::vector<float> gpu_alaw_vector;
    std::vector<float> cpu_ulaw_vector;
    std::vector<float> gpu_ulaw_vector;
    std::vector<float> cpu_g722_vector;
    std::vector<float> gpu_g722_vector;

    std::vector<float>::iterator flt_itr;

    int start_of_iterations = DEFAULT_START_OF_ITERATIONS;
    int end_of_iterations = DEFAULT_END_OF_ITERATIONS;
    int number_of_iterations = DEFAULT_NUMBER_OF_ITERATIONS;
    int number_of_averaging = DEFAULT_NUMBER_OF_AVERAGING;

    // process input arguments and generate necessary data
    if (3 == argc)
    {
        start_of_iterations = atoi(argv[1]);
        number_of_iterations = atoi(argv[2]);
    }
    else if (4==argc)
    {
        start_of_iterations = atoi(argv[1]);
        number_of_iterations = atoi(argv[2]);
        number_of_averaging = atoi(argv[3]);
    }
    else
    {
        // values are already set, print for usage for customization
        std::cout << "calculation will be done by defaut values, to customize, for example" << std::endl; 
        std::cout << "in order to calculate 13 iterations from 250 to 1024000 use..." << std::endl;
        std::cout << "--->   test 250 13" << std::endl << std::endl;
    }

    std::cout << "simulation will run with following parameters" << std::endl;
    std::cout << "start of iterations : " << start_of_iterations << "  " 
              << "end of iterations: " << end_of_iterations << "  "
              << "number of iterations : " << number_of_iterations << std::endl << std::endl;

    iteration_vector.resize(number_of_iterations);
    float step_size = (end_of_iterations - start_of_iterations) / ((float)(number_of_iterations - 1));
    int current_iteration = start_of_iterations;

    for (itr = iteration_vector.begin(); itr != iteration_vector.end(); current_iteration = current_iteration << 1)
         *itr++ = round(current_iteration);

    std::ofstream out_file;
    out_file.open("iteration_results.txt"); 

    float avg_cpu_alaw_time = 0;
    float avg_gpu_alaw_time = 0;
    float avg_cpu_ulaw_time = 0;
    float avg_gpu_ulaw_time = 0;
    float avg_cpu_g722_time = 0;
    float avg_gpu_g722_time = 0;

    cpu_alaw_vector.resize(number_of_averaging);
    gpu_alaw_vector.resize(number_of_averaging);
    cpu_ulaw_vector.resize(number_of_averaging);
    gpu_ulaw_vector.resize(number_of_averaging);
    cpu_g722_vector.resize(number_of_averaging);
    gpu_g722_vector.resize(number_of_averaging);

    unsigned int no_of_test_packet = end_of_iterations;
    unsigned int no_of_test_data = no_of_test_packet * 160;
    unsigned char* alaw_data_ptr = NULL;
    unsigned char* ulaw_data_ptr = NULL;
    unsigned char* g722_data_ptr = NULL;
    short int* cpu_decoded_alaw_data_ptr = NULL;
    short int* cpu_decoded_ulaw_data_ptr = NULL;
    short int* cpu_decoded_g722_data_ptr = NULL;	
    short int* gpu_decoded_alaw_data_ptr = NULL;
    short int* gpu_decoded_ulaw_data_ptr = NULL;
    short int* gpu_decoded_g722_data_ptr = NULL;
    unsigned int size_of_g722_band = no_of_test_packet * 45;
    int* cpu_band_ptr = NULL;
    int* gpu_band_ptr = NULL;
    int* g722_band_ptr = NULL;

    alaw_data_ptr = new unsigned char[no_of_test_data];
    ulaw_data_ptr = new unsigned char[no_of_test_data];
    g722_data_ptr = new unsigned char[no_of_test_data];	
	
    cpu_decoded_alaw_data_ptr = new short int[no_of_test_data];
    cpu_decoded_ulaw_data_ptr = new short int[no_of_test_data];
    cpu_decoded_g722_data_ptr = new short int[no_of_test_data];	
    gpu_decoded_alaw_data_ptr = new short int[no_of_test_data];
    gpu_decoded_ulaw_data_ptr = new short int[no_of_test_data];
    gpu_decoded_g722_data_ptr = new short int[no_of_test_data];
	
    cpu_band_ptr = new int[size_of_g722_band];
    gpu_band_ptr = new int[size_of_g722_band];
    g722_band_ptr = new int[size_of_g722_band];
    //std::memset(cpu_band_ptr, 0, size_of_g722_band * sizeof(int));
    //std::memset(gpu_band_ptr, 0, size_of_g722_band * sizeof(int));
	
    srand(time(NULL));

    // initialize gpu
    CudaGpuInitialize();


    std::cout << "iterations" << '\t' 
              << "g711a cpu" << '\t' << "g711a gpu" << '\t' << "speedup" << '\t' 
              << "g711u cpu" << '\t' << "g711u gpu" << '\t' << "speedup" << '\t'
              << "g722 cpu" << '\t' << "g722 gpu" << '\t' << "speedup" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------------------" << std::endl;

    out_file << "iterations" << '\t' 
             << "g711a cpu" << '\t' << "g711a gpu" << '\t' << "speedup" << '\t' 
             << "g711u cpu" << '\t' << "g711u gpu" << '\t' << "speedup" << '\t'
             << "g722 cpu" << '\t' << "g722 gpu" << '\t' << "speedup" << std::endl;

    // loop through iteration vector
    std::cout << std::fixed << std::setprecision(5);
    for (itr = iteration_vector.begin(); itr != iteration_vector.end();)
    {	
        no_of_test_packet = *itr++;
        no_of_test_data = no_of_test_packet * 160;
        size_of_g722_band = no_of_test_packet * 45;

        avg_cpu_alaw_time = 0;
        avg_gpu_alaw_time = 0;
        avg_cpu_ulaw_time = 0;
        avg_gpu_ulaw_time = 0;
        avg_cpu_g722_time = 0;
        avg_gpu_g722_time = 0;
		
        float cpu_ulaw_time = 0;
        float cpu_alaw_time = 0;
        float cpu_g722_time = 0;

        for (int k=0; k<number_of_averaging; k++)
        { 

        // generate alaw test data
        GenerateG711aData(alaw_data_ptr, no_of_test_data, COALESCED);
	
        // generate ulaw test data
        GenerateG711uData(ulaw_data_ptr, no_of_test_data, COALESCED);
	
        // generate g722 test data
        GenerateG722Data(g722_data_ptr, g722_band_ptr, no_of_test_data, COALESCED);
        std::memcpy(cpu_band_ptr, g722_band_ptr, size_of_g722_band * sizeof(int));
        std::memcpy(gpu_band_ptr, g722_band_ptr, size_of_g722_band * sizeof(int));
	
        // calculate cpu cases once since it takes too much time to complete
        //if (10  > k)
        {
        // decode ulaw using cpu
        std::clock_t cpu_ulaw_start = std::clock(); 
        G711uToPcm(ulaw_data_ptr, no_of_test_data, cpu_decoded_ulaw_data_ptr);
        std::clock_t cpu_ulaw_stop = std::clock(); 
        cpu_ulaw_time = (cpu_ulaw_stop - cpu_ulaw_start)/((float) CLOCKS_PER_SEC);
        //std::cout << cpu_ulaw_time << std::endl;

        // decode alaw using cpu
        std::clock_t cpu_alaw_start = std::clock(); 
        G711aToPcm(alaw_data_ptr, no_of_test_data, cpu_decoded_alaw_data_ptr);
        std::clock_t cpu_alaw_stop = std::clock();
        cpu_alaw_time = (cpu_alaw_stop - cpu_alaw_start)/((float) CLOCKS_PER_SEC);
        //std::cout << cpu_alaw_time << std::endl;

        // decode g722 using cpu
        std::clock_t cpu_g722_start = std::clock(); 
        G722ToPcm(g722_data_ptr, cpu_band_ptr, no_of_test_data, cpu_decoded_g722_data_ptr);
        std::clock_t cpu_g722_stop = std::clock(); 
        cpu_g722_time = (cpu_g722_stop - cpu_g722_start)/((float) CLOCKS_PER_SEC);
        }

        // decode ulaw using gpu
        std::clock_t gpu_ulaw_start = std::clock(); 
        CudaG711uToPcm(ulaw_data_ptr, no_of_test_data, gpu_decoded_ulaw_data_ptr);
        std::clock_t gpu_ulaw_stop = std::clock(); 
        float gpu_ulaw_time = (gpu_ulaw_stop - gpu_ulaw_start)/((float) CLOCKS_PER_SEC);
	
        // decode alaw using gpu
        std::clock_t gpu_alaw_start = std::clock(); 
        CudaG711aToPcm(alaw_data_ptr, no_of_test_data, gpu_decoded_alaw_data_ptr);
        std::clock_t gpu_alaw_stop = std::clock();
        float gpu_alaw_time = (gpu_alaw_stop - gpu_alaw_start)/((float) CLOCKS_PER_SEC);
	
        // decode g722 using gpu
        std::clock_t gpu_g722_start = std::clock(); 
        CudaG722ToPcm(g722_data_ptr, gpu_band_ptr, no_of_test_data, gpu_decoded_g722_data_ptr);
        std::clock_t gpu_g722_stop = std::clock(); 
        float gpu_g722_time = (gpu_g722_stop - gpu_g722_start)/((float) CLOCKS_PER_SEC);

        /*out_file  << no_of_test_packet << "\t" 
                  << cpu_alaw_time << "\t" << gpu_alaw_time << "\t"
                  << cpu_ulaw_time << "\t" << gpu_ulaw_time << "\t"
                  << cpu_g722_time << "\t" << gpu_g722_time << std::endl; */

        cpu_alaw_vector[k] = cpu_alaw_time;
        gpu_alaw_vector[k] = gpu_alaw_time;
        cpu_ulaw_vector[k] = cpu_ulaw_time;
        gpu_ulaw_vector[k] = gpu_ulaw_time;
        cpu_g722_vector[k] = cpu_g722_time;
        gpu_g722_vector[k] = gpu_g722_time;

        avg_cpu_alaw_time += (cpu_alaw_time/number_of_averaging);
        avg_gpu_alaw_time += (gpu_alaw_time/number_of_averaging);
        avg_cpu_ulaw_time += (cpu_ulaw_time/number_of_averaging);
        avg_gpu_ulaw_time += (gpu_ulaw_time/number_of_averaging);
        avg_cpu_g722_time += (cpu_g722_time/number_of_averaging);
        avg_gpu_g722_time += (gpu_g722_time/number_of_averaging);
        }

        // calculate std variance
        float std_cpu_alaw_time = 0;
        float std_gpu_alaw_time = 0;
        float std_cpu_ulaw_time = 0;
        float std_gpu_ulaw_time = 0;
        float std_cpu_g722_time = 0;
        float std_gpu_g722_time = 0;

        for (int k=0; k<number_of_averaging; k++)
        {
            std_cpu_alaw_time += pow((cpu_alaw_vector[k] - avg_cpu_alaw_time), 2.0);
            std_gpu_alaw_time += pow((gpu_alaw_vector[k] - avg_gpu_alaw_time), 2.0);
            std_cpu_ulaw_time += pow((cpu_ulaw_vector[k] - avg_cpu_ulaw_time), 2.0);
            std_gpu_ulaw_time += pow((gpu_ulaw_vector[k] - avg_gpu_ulaw_time), 2.0);
            std_cpu_g722_time += pow((cpu_g722_vector[k] - avg_cpu_g722_time), 2.0);
            std_gpu_g722_time += pow((gpu_g722_vector[k] - avg_gpu_g722_time), 2.0);       
        }

        // display results
        std::cout << no_of_test_packet << "\t\t" 
                  << avg_cpu_alaw_time << '(' << std_cpu_alaw_time << ')' << "\t" 
                  << avg_gpu_alaw_time << '(' << std_gpu_alaw_time << ')' << "\t" 
		  << avg_cpu_alaw_time/avg_gpu_alaw_time << "\t"
                  << avg_cpu_ulaw_time << '(' << std_cpu_ulaw_time << ')' << "\t" 
                  << avg_gpu_ulaw_time << '(' << std_gpu_ulaw_time << ')' << "\t" 
		  << avg_cpu_ulaw_time/avg_gpu_ulaw_time << "\t"			  
                  << avg_cpu_g722_time << '(' << std_cpu_g722_time << ')' << "\t" 
                  << avg_gpu_g722_time << '(' << std_gpu_g722_time << ')' << "\t"
				  << avg_cpu_g722_time/avg_gpu_g722_time << std::endl;
   
        out_file  << no_of_test_packet << "\t" 
                  << avg_cpu_alaw_time << '(' << std_cpu_alaw_time << ')' << "\t" 
                  << avg_gpu_alaw_time << '(' << std_gpu_alaw_time << ')' << "\t" 
		  << avg_cpu_alaw_time/avg_gpu_alaw_time << "\t"
                  << avg_cpu_ulaw_time << '(' << std_cpu_ulaw_time << ')' << "\t" 
                  << avg_gpu_ulaw_time << '(' << std_gpu_ulaw_time << ')' << "\t" 
		  << avg_cpu_ulaw_time/avg_gpu_ulaw_time << "\t"
                  << avg_cpu_g722_time << '(' << std_cpu_g722_time << ')' << "\t" 
                  << avg_gpu_g722_time << '(' << std_gpu_g722_time << ')' << "\t"
		  << avg_cpu_g722_time/avg_gpu_g722_time << std::endl; 

    }
  
    CompareArrays<short int>(cpu_decoded_alaw_data_ptr, gpu_decoded_alaw_data_ptr, no_of_test_data);
    std::cout << std::endl;

    CompareArrays<short int>(cpu_decoded_ulaw_data_ptr, gpu_decoded_ulaw_data_ptr, no_of_test_data);		
    std::cout << std::endl;

    CompareArrays<short int>(cpu_decoded_g722_data_ptr, gpu_decoded_g722_data_ptr, no_of_test_data);
    CompareArrays<int>(cpu_band_ptr, gpu_band_ptr, size_of_g722_band);
    std::cout << std::endl;

    out_file.close();
	
    if(alaw_data_ptr)
        delete []alaw_data_ptr;
    if(ulaw_data_ptr)
        delete []ulaw_data_ptr;
    if(g722_data_ptr)
        delete []g722_data_ptr;
    if(cpu_decoded_alaw_data_ptr)
        delete []cpu_decoded_alaw_data_ptr;
    if(cpu_decoded_ulaw_data_ptr)
        delete []cpu_decoded_ulaw_data_ptr;
    if(cpu_decoded_g722_data_ptr)
        delete []cpu_decoded_g722_data_ptr;
    if(cpu_band_ptr)
        delete []cpu_band_ptr;	
    if(gpu_decoded_alaw_data_ptr)
        delete []gpu_decoded_alaw_data_ptr;
    if(gpu_decoded_ulaw_data_ptr)
        delete []gpu_decoded_ulaw_data_ptr;
    if(gpu_decoded_g722_data_ptr)
        delete []gpu_decoded_g722_data_ptr;	
    if(gpu_band_ptr)
        delete []gpu_band_ptr;	
	
    return 0;
}
