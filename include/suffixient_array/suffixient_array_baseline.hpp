// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_array_baseline: Baseline suffixient array index implementation
 */

#ifndef SUFFIXIENT_ARRAY_BASELINE_HPP_
#define SUFFIXIENT_ARRAY_BASELINE_HPP_

#include <cmath>
#include <common_.hpp>
#include <sdsl/int_vector.hpp>

namespace suffixient{

template<class text_oracle>
class suffixient_array_baseline{

public:

	suffixient_array_baseline(){ this->alph = sdsl::int_vector<>(SIGMA,0,40); }

	void build(std::string input_files_basepath, text_oracle* O_, int_t len_ = 14,
	                                                        bool_t verbose = true)
	{
		{
			this->O = O_;
			this->len = len_;
		}
		std::ifstream file_alph(input_files_basepath+".mult", std::ios::binary);
		if (!file_alph.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".mult" << std::endl; exit(1);
		}
	    usafe_t a = 0, i = 0; 
	    std::vector<uchar_t> buffer(5,0);
	    while (file_alph.read(reinterpret_cast<char*>(&buffer[0]), 5))
	    {
	    	a = get_5bytes_uint(&buffer[0]);
	    	if(a > 0 and dna_to_code_table[i] > 3)
	    		{ std::cerr << "Warn: Non DNA character detected!" << std::endl; }
	    	this->alph[i] = a;
	    	i++;
	    }
	    file_alph.close();
	    {
		    usafe_t sum = 0;
		    for(i=0;i<SIGMA;++i){ usafe_t tmp = this->alph[i]; this->alph[i] = sum; sum += tmp; }
	  	}
	  	{
		  	std::ifstream file_text(input_files_basepath, std::ios::binary);
		    file_text.seekg(0, std::ios::end);
		    this->N = file_text.tellg();
		    file_text.close();
		}
	   	std::ifstream file_suff(input_files_basepath+".suff", std::ios::binary);
	   	file_suff.seekg(0, std::ios::end);
	   	this->S = file_suff.tellg()/5;
	   	file_suff.seekg(0, std::ios::beg);
	   	{
		    int_t log_n = bitsize(this->N);
		    this->Suff = sdsl::int_vector<>(S,0,log_n);
		    if(verbose)
		    	std::cout << "Suffixient array size = " << log_n << " bits per entry"
		     << std::endl;
		}
		{
			auto tmp = this->alph;
			while (file_suff.read(reinterpret_cast<char*>(&buffer[0]), 5))
			{
				a = get_5bytes_uint(&buffer[0]);
				uchar_t c = O->extract(a);
				this->Suff[tmp[c]++] = a;
			}
		} 
		if(verbose)
		{
			std::cout << "Suffixient array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S <<
			std::endl << "Suffixient array = ";
		}
	}

	usafe_t sA_size(){ return this->S; }

	int_t get_len(){ return this->len; }

	usafe_t store(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&S, sizeof(S));
		out.write((char*)&len, sizeof(len));
		w_bytes += sizeof(N) + sizeof(S) + sizeof(len);

		w_bytes += Suff.serialize(out);
		w_bytes += alph.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));
		in.read((char*)&len, sizeof(len));

		Suff.load(in);
		alph.load(in);
	}

	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix(std::string& pattern,uint_t pstart,uint_t pend) const
	{
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		int_t lcp_low, lcp_high, lcp_mid;
		plen = pend - pstart;
		low  = this->alph[pattern[pend-1]];
		high = this->alph[pattern[pend-1]+1];

		// stop if first pattern character doesn't occur in the text
		if((high - low) > 0)
			{ 
				high--;
				lcp_low = lcp_high = -1; 
				mid = (low+high)/2;
				if(plen == 1)
					return std::make_tuple(this->Suff[mid],1,false);
			}
		else
			return std::make_tuple(-1,0,true);

		while( high-low > 1 )
		{		
			auto j = O->LCS_char(pattern,pend-1,this->Suff[mid]); 
	
			if(j.first == plen)
				return std::make_tuple(this->Suff[mid],plen,false); 		
 
			if(j.second > pattern[pend-j.first-1]){
				high = mid;
				lcp_high = j.first;
			}
			else{
				low = mid;
				lcp_low = j.first;
			}
			mid = (low+high)/2;
		}

		if(lcp_low  == -1){ lcp_low  = O->LCS_char(pattern,pend-1,Suff[low]).first; }
		if(lcp_high == -1){ lcp_high = O->LCS_char(pattern,pend-1,Suff[high]).first;}

		if(lcp_low >= lcp_high){
			mid = low;
			lcp_mid = lcp_low;
		}
		else{
			mid = high;
			lcp_mid = lcp_high;
		}

		return std::make_tuple(this->Suff[mid],lcp_mid,(lcp_mid != plen));
	}

private:

	text_oracle* O;

	sdsl::int_vector<> Suff;
	sdsl::int_vector<> alph;

	uint_t S;
	usafe_t N;
	int_t len;
};

}

#endif // SUFFIXIENT_ARRAY_BASELINE_HPP_