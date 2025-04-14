#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_array_index.hpp"

template<class indexType, class oracleType>
void load_index_locate(std::string textPath, std::string oraclePath, 
                       std::string indexPath, std::string patternPath,
                       bool_t prefixArraySearch = false,
                       bool_t check_correctness = false,
                       bool_t runPrefixHeuristic = false)
{
    // initialize the index
    suffixient::suffixient_array_index
            <indexType,oracleType> suff_index;
    // construct baseline index
    std::cout << "Loading the suffixient index from " 
              << indexPath << " and " << oraclePath << std::endl;
    suff_index.load(oraclePath,indexPath);
    std::cout << "Locating the patterns in " 
              << patternPath << std::endl;

    if(check_correctness)
        suff_index.check_exact_pattern_matching_correctness(
                   textPath,patternPath,prefixArraySearch,runPrefixHeuristic);
    else
        suff_index.run_exact_pattern_matching_fasta(
                   patternPath,prefixArraySearch,runPrefixHeuristic); 
}

template<class oracleType>
void load_index_locate_ef_opt
                      (std::string textPath, std::string oraclePath, 
                       std::string indexPath, std::string patternPath,
                       bool_t check_correctness = false)
{
    // initialize the index
    suffixient::suffixient_array_index
            <suffixient::suffixient_array_elias_fano<oracleType,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>
                     ,oracleType> suff_index;
    // construct baseline index
    std::cout << "Loading the suffixient index from " 
              << indexPath << " and " << oraclePath << std::endl;
    suff_index.load(oraclePath,indexPath);
    std::cout << "Locating the patterns in " 
              << patternPath << std::endl;

    if(check_correctness)
        suff_index.check_exact_pattern_matching_correctness_ef_opt(textPath,patternPath);
    else
        suff_index.run_exact_pattern_matching_fasta_ef_opt(patternPath); 
}

void help(){

    cout << "locate [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Index files base path. Default: Empty." << endl <<
    "-t <arg>    Index variant (suffixient-array|elias-fano-opt|prefix-array). Default: suffixient-array." << endl << 
    "-o <arg>    Random access text oracle (plain-text|bitpacked-text|lz77|rlz). Default: lz77." << endl << 
    "-p <arg>    Fasta file path containing the patterns to locate. Default: Empty." << endl <<
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
        load_index_locate
        <suffixient::suffixient_array_baseline<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".sA",patternFile,false,correctness,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "plain-text")
    {
        load_index_locate
        <suffixient::suffixient_array_baseline<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".sA",patternFile,false,correctness,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "bitpacked-text")
    {
        load_index_locate
        <suffixient::suffixient_array_baseline<suffixient::bitpacked_text_oracle>,
         suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".sA",patternFile,false,correctness,true);
    }
    else if(indexType == "suffixient-array" and oracleType == "rlz")
    {
        load_index_locate
        <suffixient::suffixient_array_baseline<RLZ_DNA<>>,
         RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".sA",patternFile,false,correctness,true);
    }
    else if(indexType == "prefix-array" and oracleType == "lz77")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".pa",patternFile,true,correctness,false);
    }
    else if(indexType == "prefix-array" and oracleType == "plain-text")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".pa",patternFile,true,correctness,false);
    }
    else if(indexType == "prefix-array" and oracleType == "bitpacked-text")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<suffixient::bitpacked_text_oracle>,
         suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".pa",patternFile,true,correctness,false);
    }
    else if(indexType == "prefix-array" and oracleType == "rlz")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<RLZ_DNA<>>,
         RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".pa",patternFile,true,correctness,false);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "lz77")
    {
        load_index_locate_ef_opt<lz77::LZ77_compressed_text>
        (inputPath,inputPath+".lz77",inputPath+".opt_sA",patternFile,correctness);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "plain-text")
    {
        load_index_locate_ef_opt<suffixient::uncompressed_text_oracle>
        (inputPath,inputPath,inputPath+".opt_sA",patternFile,correctness);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "bitpacked-text")
    {
        load_index_locate_ef_opt<suffixient::bitpacked_text_oracle>
        (inputPath,inputPath,inputPath+".opt_sA",patternFile,correctness);
    }
    else if(indexType == "elias-fano-opt" and oracleType == "rlz")
    {
        load_index_locate_ef_opt<RLZ_DNA<>>
        (inputPath,inputPath+".rlz",inputPath+".opt_sA",patternFile,correctness);
    }
    else{
        std::cout << "Not a valid input parameter configuration..."  <<
                     " Please use the -h flag for more information." << std::endl;
        exit(1);
    }

    return 0;
}