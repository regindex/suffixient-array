#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

#include "suffixient_DNA_index.hpp"

void help(){

    cout << "build_DNA_index [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-i <arg>    Input files base path. Default: Empty." << endl <<
    "-t <arg>    Index type (baseline|elias-fano|heuristic). Default: baseline." << endl << 
    "-l <arg>    Maximum string length to store in the elias-fano data structure. Default: 14." << endl << 
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

    std::string inputPath, indexType = "baseline";
    int_t len = 14;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hvi:t:l:")) != -1)
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
            case 'l':
                len = std::atoi(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    if(indexType == "baseline")
    {
        // initialize the index
        suffixient::suffixient_DNA_index
                <suffixient::suffixient_array_baseline<lz77::LZ77_Kreft_Navarro_index>,
                lz77::LZ77_Kreft_Navarro_index> DNA_index;
        // construct baseline index
        DNA_index.build(inputPath, len);
    }
    else{
        std::cout << "Not yet implemented..." << std::endl;
        exit(1);
    }

    return 0;
}