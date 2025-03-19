#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <limits>
#include <unistd.h>

using namespace std;

#include "suffixient_zfast_index.hpp"

void help(){

    cout << "build_test_index [options]" << endl <<
    "Options:" << endl <<
    "-h          Print usage info." << endl <<
    "-t <arg>    Input text file path. Default: Empty." << endl << 
    "-s <arg>    Suffixient set file path. Default: Empty." << endl << 
    "-l <arg>    Supermaximal string lengths file path. Default: Empty." << endl << 
    "-p <arg>    Fasta file path containing the sequences for MEMs searching procedure. Default: Empty." << endl << 
    //"-o <arg>    Store output to file using 64-bits unsigned integers. If not specified, output is streamed to standard output in human-readable format." << endl <<
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

    string outputFile, textFile, suffixientFile, lengthsFile, patternsFile;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "hvp:l:s:t:o:")) != -1)
    {
        switch (opt){
            case 'h':
                help();
            break;
            case 'v':
                verbose=true;
            break;
            case 'o':
                outputFile = string(optarg);
            break;
            case 't':
                textFile = string(optarg);
            break;
            case 's':
                suffixientFile = string(optarg);
            break;
            case 'l':
                lengthsFile = string(optarg);
            break;
            case 'p':
                patternsFile = string(optarg);
            break;
            default:
                help();
            return -1;
        }
    }

    // initialize the index
    suffixient::suffixient_index
            <ctriepp::CTriePP<int_t>,lz77::LZ77_Kreft_Navarro_index> index{};
    // build the index
    std::pair<safe_t,safe_t> res;
    if(lengthsFile.size() == 0)
        res = index.build(textFile,suffixientFile,verbose);
    else
        res = index.build(textFile,suffixientFile,lengthsFile,verbose);

    std::ofstream output(textFile+".stats");

    output << "Index construction completed for " << textFile << ":" << std::endl
              << "Total inserted keywords = " << res.second << std::endl
              << "Total inserted characters = " << res.first << std::endl;

    output.close();

    if(patternsFile.size() > 0)
        index.locate_fasta(patternsFile);

    return 0;
}