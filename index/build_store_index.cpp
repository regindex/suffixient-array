#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_sA_index.hpp"

int_t compute_stored_strings_length(std::string inputPath, uint_t addSpace)
{
    std::ifstream file_text(inputPath, std::ios::binary);
    file_text.seekg(0, std::ios::end);
    usafe_t N = file_text.tellg();
    file_text.seekg(0, std::ios::beg);
    file_text.close();

    std::ifstream file_suff(inputPath+".suff", std::ios::binary);
    file_suff.seekg(0, std::ios::end);
    usafe_t S = file_suff.tellg()/5;
    file_suff.seekg(0, std::ios::beg);
    file_suff.close();

    usafe_t sA_size = S * ceil(log2(N));
    double additional = ((double)sA_size/100) * addSpace;

    usafe_t bv_size = ((S+1)/64 + 1) * 64; bv_size += bv_size * 0.2;
    while(bv_size > additional)
    {
        addSpace += 1;
        additional = ((double)sA_size/100) * addSpace;
    }

    int_t len = 1;
    while(true)
    {
        usafe_t S_ = std::min(static_cast<usafe_t>(pow(4,len)),static_cast<usafe_t>(S));
        usafe_t ef_size = S_ * (2 + log2(pow(4,len)/S_));
        if((bv_size + ef_size + 128) > additional){ len--; break; }
        len++;
    }
    if(len == 0) len++;
    
    return len;
}

template<class indexType, class oracleType>
void construct_store_index(std::string inputPath, std::string indexPath,
                                                     int_t addSpace = 0)
{
    /*
    int_t storedLen = 0;
    if(addSpace > 0){
        storedLen = compute_stored_strings_length(inputPath, addSpace);
        std::cout << "Storing substring of length " << storedLen << " preceding"
                  << " the suffixient positions" << std::endl; 
    }
    else{ storedLen = 10; }
    */
    int_t storedLen = 14;

    std::cout << "Contructing and storing the suffixient index to " 
              << indexPath << std::endl;
    // initialize the index
    suffixient::suffixient_sA_index
            <indexType,oracleType> suff_index;
    // construct baseline index
    suff_index.build(inputPath,storedLen);
    std::cout << "Index size = " << 
    suff_index.store(indexPath)
    << " bytes" << std::endl;
}

void help(){

    cout << "build_DNA_index [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Input files base path. Default: Empty." << endl <<
    "-t <arg>    Index type (baseline|elias-fano|prefix-array). Default: baseline." << endl << 
    "-o <arg>    Text oracle (plain|lz77|Hk). Default: lz77." << endl << 
    "-l <arg>    Maximum additional space (percentage) on top of the suffixient array. Default: 10%." << endl << 
    "-v          Activate verbosity mode. Default: false." << endl;
    exit(0);
} 

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cerr << "Wrong number of parameters... See the help messagge:" << endl;
        help();
        exit(1);
    }

    std::string inputPath, indexType = "baseline", oracleType = "lz77";
    int_t additional_space = 10;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hvo:i:t:l:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'v':
                verbose = true;
            break;
            case 'i':
                inputPath = string(optarg);
            break;
            case 't':
                indexType = string(optarg);
            break;
            case 'o':
                oracleType = string(optarg);
            break;
            case 'l':
                additional_space = std::atoi(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(indexType == "baseline" and oracleType == "lz77")
    {
        construct_store_index
        <suffixient::suffixient_array_baseline<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".bai");
    }
    else if(indexType == "baseline" and oracleType == "plain")
    {
        construct_store_index
        <suffixient::suffixient_array_baseline<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".bai");
    }
    else if(indexType == "baseline" and oracleType == "Hk")
    {
        construct_store_index
        <suffixient::suffixient_array_baseline<Hk_Ferragina_Venturini<>>,
         Hk_Ferragina_Venturini<>>
        (inputPath,inputPath+".bai");
    }
    else if(indexType == "prefix-array" and oracleType == "lz77")
    {
        construct_store_index
        <suffixient::suffix_array_binary_search<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".pai");
    }
    else if(indexType == "prefix-array" and oracleType == "plain")
    {
        construct_store_index
        <suffixient::suffix_array_binary_search<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".pai");
    }
    else if(indexType == "prefix-array" and oracleType == "Hk")
    {
        construct_store_index
        <suffixient::suffix_array_binary_search<Hk_Ferragina_Venturini<>>,
         Hk_Ferragina_Venturini<>>
        (inputPath,inputPath+".pai");
    }
    else if(indexType == "elias-fano" and oracleType == "lz77")
    {
        //std::string extension = "." + std::to_string(additional_space) + "efi";
        std::string extension = ".efi";
        construct_store_index
        <suffixient::suffixient_array_elias_fano<lz77::LZ77_compressed_text,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     lz77::LZ77_compressed_text>
        (inputPath,inputPath+extension,additional_space);
    }
    else if(indexType == "elias-fano" and oracleType == "plain")
    {
        //std::string extension = "." + std::to_string(additional_space) + "efi";
        std::string extension = ".efi";
        construct_store_index
        <suffixient::suffixient_array_elias_fano<suffixient::uncompressed_text_oracle,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+extension,additional_space);
    }
    else if(indexType == "elias-fano" and oracleType == "Hk")
    {
        //std::string extension = "." + std::to_string(additional_space) + "efi";
        std::string extension = ".efi";
        construct_store_index
        <suffixient::suffixient_array_elias_fano<Hk_Ferragina_Venturini<>,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     Hk_Ferragina_Venturini<>>
        (inputPath,inputPath+extension,additional_space);
    }
    else{
        std::cout << "Not yet implemented..." << std::endl;
        exit(1);
    }

    return 0;
}