// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_array_baseline: DESCRIPTION
 */

#ifndef SUFFIXIENT_ARRAY_BASELINE_HPP_
#define SUFFIXIENT_ARRAY_BASELINE_HPP_

#include <common.hpp>
#include <sdsl/int_vector.hpp>

namespace suffixient{

template<class text_oracle>
class suffixient_array_baseline{

public:

	suffixient_array_baseline(){ this->alph = std::vector<usafe_t>(5,0); }

	void build(std::string input_files_basepath, const uint8_t len_)
	{
		this->len = len_;
		if(this->len > 32)
		{ 
			std::cerr << "Length must be smaller or equal than 32!" << std::endl; 
			exit(1);
		}
		std::ifstream file_alph(input_files_basepath+".alph", std::ios::binary);
		if (!file_alph.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".alph" << std::endl; exit(1);
		}
	    usafe_t a, i = 0;
	    while (file_alph.read(reinterpret_cast<char*>(&a), 5))
	    {
	    	if(a > 0 and dna_to_code_table[i] > 3)
	    		{ std::cerr << "Error: Non DNA character detected!" << std::endl; 
	    		  //exit(1); }
	    		}
	    	else{ this->alph[dna_to_code_table[i]] = a; }
	    	i++;
	    }
	    file_alph.close();
	    {
		    usafe_t sum = 0;
		    for(i=0;i<5;++i){ usafe_t tmp = alph[i]; this->alph[i] = sum; sum += tmp; }
	  	}
	  	/*
	    // Open the input files in binary mode
	  	std::ifstream fileT(input_files_basepath, std::ios::binary);
	    std::ifstream fileS(input_files_basepath+".suff", std::ios::binary);
	    std::ifstream fileL(input_files_basepath+".lcs", std::ios::binary);
	    if (!fileT.is_open() or !fileS.is_open() or !fileL.is_open()) {
	        std::cerr << "Error: Could not open input files... " 
	        					<< input_files_basepath         << ", " 
	        					<< input_files_basepath+".suff" << ", " 
	        					<< input_files_basepath+".lcs" << std::endl; exit(1); }

	    Ulong j = 0;
	    while (fileS.read(reinterpret_cast<char*>(&a), 5) and
	    	   fileL.read(reinterpret_cast<char*>(&i), 5))
	    {
	    	i += 1; std::cout << "i: " << i << std::endl;
	    	Ulong start = a-i+1;
	    	std::string buffer(i,0);
	    	{
	    		fileT.clear();
	    		fileT.seekg(start, std::ios::beg);
	    		fileT.read(&buffer[0], i);
	    	}
	    	//for(int l=0;l<buffer.size();++l)
	    	//	std::cout << char(buffer[buffer.size()-l-1]) << std::endl;

	    	std::cout << buffer << std::endl;
	    	Long pi = 0;
	    	set_uint(reinterpret_cast<uint8_t*>(&pi),buffer);
	    	std::cout << "-->" << pi << std::endl;

	    	if(j++ > 2) exit(1);
	    }
	    
	  	uint32_t pi = 0; 
	  	std::string buff = "AAA";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "AAT";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "AT";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "GA";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "GT";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "TAA";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "TAG";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;
	  	buff = "TAT";
	  	set_uint(reinterpret_cast<uint8_t*>(&pi),buff);
	  	std::cout << "--> " << pi << std::endl;
	  	pi = 0;

	  	// if(str1 ^ str2 != 0)
	  	//	int LCP = __builtin_clz(str1 ^ str2) // number of leading zeros after xor
	    exit(1);
	    */
	}

/*
16 AAA  00000000000000000000000000000000    0
15 AAT  00001100000000000000000000000000    201326592
11 AT   00110000000000000000000000000000    805306368
17 GA   10000000000000000000000000000000    2147483648
 8 GT   10110000000000000000000000000000    2952790016  
 5 TAA  11000000000000000000000000000000    3221225472
10 TAG  11001000000000000000000000000000    3355443200
 7 TAT  11001100000000000000000000000000    3422552064
*/

private:
	//
	text_oracle* O;
	//
	sdsl::int_vector<> S;
	//
	std::vector<usafe_t> alph;
	//
	int_t len;
};

}

#endif // SUFFIXIENT_ARRAY_BASELINE_HPP_