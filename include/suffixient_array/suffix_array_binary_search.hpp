// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_array_binary_search: index performing binary search on the full prefix array
 */

#ifndef SUFFIX_ARRAY_BINARY_SEARCH_HPP_
#define SUFFIX_ARRAY_BINARY_SEARCH_HPP_

#include <cmath>
#include <common.hpp>
#include <sdsl/int_vector.hpp>

namespace suffixient{

template<class text_oracle>
class suffix_array_binary_search{

public:

	suffix_array_binary_search(){ this->alph = sdsl::int_vector<>(SIGMA,0,40); }

	int_t get_len(){ return this->len; }

	void build(std::string input_files_basepath, text_oracle* O_, int_t len_ = 14,
	                                                       bool_t verbose = true)
	{
		{
			this->O = O_;
			this->len = len_;
		}
	  	{
		  	std::ifstream file_text(input_files_basepath, std::ios::binary);
		    file_text.seekg(0, std::ios::end);
		    this->N = this->S = file_text.tellg();
		    file_text.seekg(0, std::ios::beg);
		    uchar_t c;
			while (file_text.read(reinterpret_cast<char*>(&c), 1)){ this->alph[c]++; }
		    file_text.close();
		}
		usafe_t a = 0, i = 0;
	    {
		    usafe_t sum = 0;
		    for(i=0;i<SIGMA;++i){ usafe_t tmp = this->alph[i]; this->alph[i] = sum; sum += tmp; }
	  	}
	   	{
		    int_t log_n = bitsize(this->S);
		    this->SA = sdsl::int_vector<>(this->S,0,log_n);
		    if(verbose)
		    	std::cout << "Suffix array size = " << log_n << " bits per entry"
		     << std::endl;
		}
	   	std::ifstream file_suff(input_files_basepath+".pa", std::ios::binary);
		if (!file_suff.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".pa" << std::endl; exit(1);
		}
		{
			i = 0;
			auto tmp = this->alph;
			std::vector<uchar_t> buffer(5,0);
			while (file_suff.read(reinterpret_cast<char*>(&buffer[0]), 5))
			{
				a = get_5bytes_uint(&buffer[0])-1;
				this->SA[i++] = a;
			}
		}
		file_suff.close();
		if(verbose)
		{
			std::cout << "Suffixient array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S <<
			std::endl << "Suffixient array = ";
		}
	}

	usafe_t store(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&S, sizeof(S));
		w_bytes += sizeof(N) + sizeof(S);

		w_bytes += SA.serialize(out);
		w_bytes += alph.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&S, sizeof(S));

		SA.load(in);
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
					return std::make_tuple(this->SA[mid],1,false);
			}
		else
			return std::make_tuple(-1,0,true);

		while( high-low > 1 )
		{		
			auto j = O->LCS_char(pattern,pend-1,this->SA[mid]); 
	
			if(j.first == plen){ return std::make_tuple(this->SA[mid],plen,false); }

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

		if(lcp_low  == -1){ lcp_low  = O->LCS_char(pattern,pend-1,SA[low]).first; }
		if(lcp_high == -1){ lcp_high = O->LCS_char(pattern,pend-1,SA[high]).first;}

		if(lcp_low >= lcp_high){
			mid = low;
			lcp_mid = lcp_low;
		}
		else{
			mid = high;
			lcp_mid = lcp_high;
		}

		return std::make_tuple(this->SA[mid],lcp_mid,(lcp_mid != plen));
	}

private:

	text_oracle* O;

	sdsl::int_vector<> SA;
	sdsl::int_vector<> alph;

	uint_t S;
	usafe_t N;
	int_t len;
};

}

#endif