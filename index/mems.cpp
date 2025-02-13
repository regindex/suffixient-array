#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_sA_index.hpp"

template<class indexType, class oracleType>
void load_index_find_mems(std::string textPath, 
                          std::string indexPath, std::string patternPath,
                          bool_t prefixArraySearch = false)
{
    // initialize the index
    suffixient::suffixient_sA_index
            <indexType,oracleType> suff_index;
    // construct baseline index
    suff_index.load(textPath,indexPath);
    std::cout << "Loading the suffixient index from " 
              << indexPath << std::endl;

    suff_index.locate_MEMs_fasta(patternPath);
}

void help(){

    cout << "locate [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Index files base path. Default: Empty." << endl <<
    "-t <arg>    Index type (baseline|elias-fano|heuristic|sa). Default: baseline." << endl << 
    "-o <arg>    Text oracle (uncompressed|lz77). Default: lz77." << endl << 
    "-p <arg>    Fasta file path containing the patterns to locate. Default: Empty." << endl <<
    //"-c          Check output correctness. Default: False." << endl << 
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

    std::string inputPath, patternFile, indexType = "baseline", oracleType = "lz77";
    bool verbose = false, correctness = false;

    int opt;
    while ((opt = getopt(argc, argv, "hvo:i:t:p:c")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'v':
                verbose = true;
            break;
            case 'c':
                correctness = true;
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
            case 'p':
                patternFile = string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(indexType == "baseline" and oracleType == "lz77")
    {
        load_index_find_mems
        <suffixient::suffixient_array_baseline<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".basi",patternFile,false);
    }
    else if(indexType == "baseline" and oracleType == "uncompressed")
    {
        load_index_find_mems
        <suffixient::suffixient_array_baseline<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".basi",patternFile,false);
    }
    else if(indexType == "sa" and oracleType == "lz77")
    {
        load_index_find_mems
        <suffixient::suffix_array_binary_search<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".sai",patternFile,true);
    }
    else if(indexType == "sa" and oracleType == "uncompressed")
    {
        load_index_find_mems
        <suffixient::suffix_array_binary_search<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".sai",patternFile,true);
    }
    else if(indexType == "elias-fano" and oracleType == "lz77")
    {
        load_index_find_mems
        <suffixient::suffixient_array_elias_fano<lz77::LZ77_compressed_text,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".efi",patternFile,false);
    }
    else if(indexType == "elias-fano" and oracleType == "uncompressed")
    {
        load_index_find_mems
        <suffixient::suffixient_array_elias_fano<suffixient::uncompressed_text_oracle,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".efi",patternFile,false);
    }
    else{
        std::cout << "Not yet implemented..." << std::endl;
        exit(1);
    }

    return 0;
}