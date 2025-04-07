#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_array_index.hpp"

template<class indexType, class oracleType>
void load_index_mems(std::string textPath, std::string oraclePath, 
                       std::string indexPath, std::string patternPath,
                       bool_t runPrefixHeuristic = false,
                       bool_t runEFopt = false)
{
    // initialize the index
    suffixient::suffixient_array_index
            <indexType,oracleType> suff_index;
    // construct baseline index
    std::cout << "Loading the suffixient index from " 
              << indexPath << " and " << oraclePath << std::endl;
    suff_index.load(oraclePath,indexPath);
    std::cout << "Finding MEMs for the patterns in " 
              << patternPath << std::endl;

    //if(check_correctness)
    //    suff_index.check_exact_pattern_matching_correctness(
    //               textPath,patternPath,prefixArraySearch,runPrefixHeuristic);
    //else
    suff_index.run_MEMs_finding_fasta(
                   patternPath,runPrefixHeuristic,runEFopt); 
}

void help(){

    cout << "locate [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Index files base path. Default: Empty." << endl <<
    "-t <arg>    Index variant (suffixient-array|elias-fano-opt|prefix-array). Default: suffixient-array." << endl << 
    "-o <arg>    Random access text oracle (plain-text|bitpacked-text|lz77|rlz). Default: lz77." << endl << 
    "-p <arg>    Fasta file path containing the patterns. Default: Empty." << endl <<
    "-c          Check output correctness. Default: False." << endl << 
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

    std::string inputPath, patternFile, indexType = "suffixient-array", oracleType = "lz77";
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

    if(indexType == "suffixient-array" and oracleType == "lz77")
    {
        load_index_mems
        <suffixient::suffixient_array_baseline<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".bai",patternFile,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "plain-text")
    {
        load_index_mems
        <suffixient::suffixient_array_baseline<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".bai",patternFile,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "bitpacked-text")
    {
        load_index_mems
        <suffixient::suffixient_array_baseline<suffixient::bitpacked_text_oracle>,
         suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".bai",patternFile,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "rlz")
    {
        load_index_mems
        <suffixient::suffixient_array_baseline<RLZ_DNA<>>,
         RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".bai",patternFile,true);
    }
    else if(indexType == "prefix-array" and oracleType == "lz77")
    {
        load_index_mems
        <suffixient::suffix_array_binary_search<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".pai",patternFile,true);
    }
    else if(indexType == "prefix-array" and oracleType == "plain-text")
    {
        load_index_mems
        <suffixient::suffix_array_binary_search<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".pai",patternFile,true);
    }
    else if(indexType == "prefix-array" and oracleType == "bitpacked-text")
    {
        load_index_mems
        <suffixient::suffix_array_binary_search<suffixient::bitpacked_text_oracle>,
         suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".pai",patternFile,true);
    }
    else if(indexType == "prefix-array" and oracleType == "rlz")
    {
        load_index_mems
        <suffixient::suffix_array_binary_search<RLZ_DNA<>>,
         RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".pai",patternFile,true);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "lz77")
    {
        load_index_mems
        <suffixient::suffixient_array_elias_fano<lz77::LZ77_compressed_text,
                                                 suffixient::elias_fano_bitvector,
                                                 suffixient::succinct_bitvector>,
         lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".efi",patternFile,true,true);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "plain-text")
    {
        load_index_mems
        <suffixient::suffixient_array_elias_fano<suffixient::uncompressed_text_oracle,
                                                 suffixient::elias_fano_bitvector,
                                                 suffixient::succinct_bitvector>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".efi",patternFile,true,true);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "bitpacked-text")
    {
        load_index_mems
        <suffixient::suffixient_array_elias_fano<suffixient::bitpacked_text_oracle,
                                                 suffixient::elias_fano_bitvector,
                                                 suffixient::succinct_bitvector>,
         suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".efi",patternFile,true,true);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "rlz")
    {
        load_index_mems
        <suffixient::suffixient_array_elias_fano<RLZ_DNA<>,
                                                 suffixient::elias_fano_bitvector,
                                                 suffixient::succinct_bitvector>,
         RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".efi",patternFile,true,true);
    }
    else{
        std::cout << "Not yet implemented..." << std::endl;
        exit(1);
    }
    
    return 0;
}