#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_sA_index.hpp"

template<class indexType, class oracleType>
void load_index_locate(std::string textPath, 
                       std::string indexPath, std::string patternPath,
                       bool_t prefixArraySearch = false,
                       bool_t check_correctness = false,
                       bool_t runPrefixHeuristic = false)
{
    // initialize the index
    suffixient::suffixient_sA_index
            <indexType,oracleType> suff_index;
    // construct baseline index
    std::cout << "Loading the suffixient index from " 
              << indexPath << std::endl;
    suff_index.load(textPath,indexPath);
    std::cout << "Locating the patterns in " 
              << indexPath << std::endl;

    if(check_correctness)
        suff_index.check_exact_pattern_matching_correctness(patternPath,prefixArraySearch,runPrefixHeuristic);
    else
        suff_index.run_exact_pattern_matching_fasta(patternPath,prefixArraySearch,runPrefixHeuristic); 
}

void help(){

    cout << "locate [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Index files base path. Default: Empty." << endl <<
    "-t <arg>    Index type (baseline|elias-fano|prefix-array). Default: baseline." << endl << 
    "-o <arg>    Text oracle (plain|lz77|Hk). Default: lz77." << endl << 
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
        load_index_locate
        <suffixient::suffixient_array_baseline<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".bai",patternFile,false,correctness,true);
    }
    else if(indexType == "baseline" and oracleType == "plain")
    {
        load_index_locate
        <suffixient::suffixient_array_baseline<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".bai",patternFile,false,correctness,true);
    }
    else if(indexType == "baseline" and oracleType == "Hk")
    {
        load_index_locate
        <suffixient::suffixient_array_baseline<Hk_Ferragina_Venturini<>>,
         Hk_Ferragina_Venturini<>>
        (inputPath+".hkfv",inputPath+".bai",patternFile,false,correctness,true);
    }
    else if(indexType == "prefix-array" and oracleType == "lz77")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<lz77::LZ77_compressed_text>,
               lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".pai",patternFile,true,correctness);
    }
    else if(indexType == "prefix-array" and oracleType == "plain")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<suffixient::uncompressed_text_oracle>,
         suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".pai",patternFile,true,correctness);
    }
    else if(indexType == "prefix-array" and oracleType == "Hk")
    {
        load_index_locate
        <suffixient::suffix_array_binary_search<Hk_Ferragina_Venturini<>>,
         Hk_Ferragina_Venturini<>>
        (inputPath+".hkfv",inputPath+".pai",patternFile,true,correctness);
    }
    else if(indexType == "elias-fano" and oracleType == "lz77")
    {
        load_index_locate
        <suffixient::suffixient_array_elias_fano<lz77::LZ77_compressed_text,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     lz77::LZ77_compressed_text>
        (inputPath+".lz77",inputPath+".efi",patternFile,false,correctness,true);
    }
    else if(indexType == "elias-fano" and oracleType == "plain")
    {
        load_index_locate
        <suffixient::suffixient_array_elias_fano<suffixient::uncompressed_text_oracle,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     suffixient::uncompressed_text_oracle>
        (inputPath,inputPath+".efi",patternFile,false,correctness,true);
    }
    else if(indexType == "elias-fano" and oracleType == "Hk")
    {
        load_index_locate
        <suffixient::suffixient_array_elias_fano<Hk_Ferragina_Venturini<>,
                     suffixient::elias_fano_bitvector,suffixient::succinct_bitvector>,
                     Hk_Ferragina_Venturini<>>
        (inputPath+".hkfv",inputPath+".efi",patternFile,false,correctness,true);
    }
    else{
        std::cout << "Not yet implemented..." << std::endl;
        exit(1);
    }
    // .sdi

    return 0;
}