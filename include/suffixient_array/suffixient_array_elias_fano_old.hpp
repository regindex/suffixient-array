// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_array_baseline: DESCRIPTION
 */

#ifndef SUFFIXIENT_ARRAY_ELIAS_FANO_HPP_
#define SUFFIXIENT_ARRAY_ELIAS_FANO_HPP_

#include <cmath>
#include <common.hpp>
#include <sdsl/int_vector.hpp>

#include <elias_fano_bitvector.hpp>
#include <succinct_bitvector.hpp>

namespace suffixient{

template<class text_oracle, class elias_fano_ds, class bitvector>
class suffixient_array_elias_fano{

public:

	suffixient_array_elias_fano(){}

	void build(std::string input_files_basepath, text_oracle* O_, 
		                             int_t len_ = 1, bool_t verbose = true)
	{
		usafe_t a = 0, i = 0, l = 0; 
	    std::vector<uchar_t> buffer(STORE_SIZE,0);
		{
			O = O_;
		}
		this->len = len_;
		if(this->len > 16 and sizeof(uint_t) == 4)
		{ 
			std::cerr << "The supermaximal substrings length must be smaller or equal than 16 when using words of 4 bytes!" << std::endl; 
			std::cerr << "Please, run the 64 bits executable..." << std::endl; 
			exit(1);
		}
		if(this->len > 32)
		{ 
			std::cerr << "The supermaximal substrings length must be smaller or equal than 32!" << std::endl; 
			exit(1);
		}
		std::vector<usafe_t> alph = std::vector<usafe_t>(SIGMA_DNA+1,0);
		std::ifstream file_alph(input_files_basepath+".mult", std::ios::binary);
		if (!file_alph.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".mult" << std::endl; exit(1);
		}
	    while (file_alph.read(reinterpret_cast<char*>(&buffer[0]), STORE_SIZE))
	    {
	    	a = get_5bytes_uint(&buffer[0]);
	    	if(a > 0 and dna_to_code_table[i] > 3)
	    		{ std::cerr << "Error: Non DNA character detected!" << std::endl; 
	    		  exit(1);
	    		}
	    	alph[dna_to_code_table[i]] = a;
	    	i++;
	    }
	    file_alph.close();
	    {
		    usafe_t sum = 0;
		    for(i=0;i<SIGMA_DNA+1;++i)
		    	{ usafe_t tmp = alph[i]; alph[i] = sum; sum += tmp; }
	  	}
	  	{
		  	std::ifstream file_text(input_files_basepath, std::ios::binary);
		    file_text.seekg(0, std::ios::end);
		    this->N = file_text.tellg();
		    file_text.close();
		}
	   	std::ifstream file_suff(input_files_basepath+".suff", std::ios::binary);
		if (!file_suff.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".suff" << std::endl; exit(1);
		}
	   	file_suff.seekg(0, std::ios::end);
	   	usafe_t S = file_suff.tellg()/5;
	   	file_suff.seekg(0, std::ios::beg);
	   	sdsl::int_vector<> LCS;
	   	{
		    int_t log_n = bitsize(this->N);
		    this->Suff = sdsl::int_vector<>(S,0,log_n);
		    LCS = sdsl::int_vector<>(S,0,log_n);
		    if(verbose)
		    	std::cout << "Suffixient array size = " << log_n << " bits per entry"
		     << std::endl;
		}
		std::ifstream file_lcs(input_files_basepath+".lcs", std::ios::binary);
		if (!file_lcs.is_open())
		{
	    std::cerr << "Error: Could not open " << 
	    			 			 input_files_basepath+".lcs" << std::endl; exit(1);
		}
		{
			while (file_suff.read(reinterpret_cast<char*>(&buffer[0]), STORE_SIZE))
			{
				a = get_5bytes_uint(&buffer[0]);
	    		file_lcs.read(reinterpret_cast<char*>(&buffer[0]), STORE_SIZE);
	    		l = get_5bytes_uint(&buffer[0]);

				uchar_t c = O->extract(a);
				this->Suff[alph[dna_to_code_table[c]]] = a;
				LCS[alph[dna_to_code_table[c]]++] = l;
			}
		} 
		if(verbose)
		{
			std::cout << "Suffixient array size = " << this->Suff.size() <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/this->Suff.size() <<
			std::endl << "Suffixient array = " << std::endl;
			for(i=0;i<this->Suff.size();++i)
				std::cout << this->Suff[i] << " " << LCS[i] << std::endl;
			std::cout << std::endl;
		}
		{
			std::vector<uint_t> onset;
			uint_t offset = 0;
			std::ifstream file_text(input_files_basepath, std::ios::binary);
			if (!file_text.is_open())
			{
		    std::cerr << "Error: Could not open " << 
		    			 			 input_files_basepath << std::endl; exit(1);
			}
			safe_t prev = -1;
			std::vector<bool_t> tmp_bv(this->Suff.size()+1,0);
			for(i=0;i<S;++i)
			{
				if(LCS[i] < this->len)
				{
					tmp_bv[i] = 1; 
					if(prev != -1)
					{
						std::string text_buffer(this->len,'A');
						safe_t beg = std::max(static_cast<safe_t>(0),prev-this->len+1);
						safe_t len_ = std::min(static_cast<safe_t>(this->len),prev+1);

						file_text.seekg(beg, std::ios::beg);
						file_text.read(&text_buffer[this->len-len_], len_);
					  	file_text.clear();

					  	//std::cout << text_buffer << " ";
					  	set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&offset),text_buffer);
					  	//std::cout << offset << std::endl;
						onset.push_back(offset);
						offset = 0;
					}
					prev = this->Suff[i];
				}
			}
			tmp_bv[this->Suff.size()] = 1;

			std::string text_buffer(this->len,'A');
			safe_t beg = std::max(static_cast<safe_t>(0),prev-this->len+1);
			safe_t len_ = std::min(static_cast<safe_t>(this->len),prev+1);
			file_text.seekg(beg, std::ios::beg);
			file_text.read(&text_buffer[this->len-len_], len_);
			//std::cout << text_buffer << " ";
			set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&offset),text_buffer);
			//std::cout << offset << std::endl;
			onset.push_back(offset);

			file_text.close();

			this->bv = bitvector(tmp_bv);
			/*
			std::cout << "ef size = " << pow(SIGMA_DNA,this->len) << std::endl;
			for(usafe_t k=0;k<onset.size();++k)
				std::cout << onset[k] << std::endl; */
			this->ef = elias_fano_ds(onset,pow(SIGMA_DNA,this->len));
			/*
			for(usafe_t k=0;k<this->Suff.size()+1;++k)
				std::cout << bv[k];
			std::cout << std::endl;*/
			file_text.close();
		}

	    file_suff.close();
	    file_lcs.close();
	}

		/*
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
	    		{ std::cerr << "Warn: Non DNA character detected!" << std::endl; 
	    		  exit(1);
	    		}
	    	else{ this->alph[dna_to_code_table[i]] = a; }
	    	i++;
	    }
	    file_alph.close();
	    {
		    usafe_t sum = 0;
		    for(i=0;i<5;++i){ usafe_t tmp = alph[i]; this->alph[i] = sum; sum += tmp; }
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
				this->Suff[tmp[dna_to_code_table[c]]++] = a;
			}
		}
		if(verbose)
		{
			std::cout << "Suffixient array size = " << this->S <<
			std::endl << "Text size = " << this->N <<
			std::endl << "N/S = " << double(this->N)/S <<
			std::endl << "Suffixient array = ";
			for(i=0;i<this->S;++i)
				std::cout << this->Suff[i] << " ";
			std::cout << std::endl;
		}
		*/

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

	usafe_t sA_size(){ return this->S; }

	int_t get_len(){ return this->len; }

	usafe_t store(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		out.write((char*)&N, sizeof(N));
		out.write((char*)&len, sizeof(len));
		w_bytes += sizeof(N) + sizeof(len);

		w_bytes += Suff.serialize(out);
		w_bytes += ef.serialize(out);
		w_bytes += bv.serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* O_)
	{
		O = O_;
		in.read((char*)&N, sizeof(N));
		in.read((char*)&len, sizeof(len));

		Suff.load(in);
		ef.load(in);
		bv.load(in);
	}

	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix(std::string& pattern,uint_t pstart,uint_t pend) const
	{
		//pattern[pstart] = 'C';
		//pattern[pstart+1] = 'G';
		//pend++;
		//pend++;
		//pend++;

		std::cout << "this->len = " << this->len << std::endl;
		std::cout << "LOCate_longest_prefix = " << pattern.substr(pstart,pend-pstart) << std::endl;

		uint_t search = 0;
		uint_t plen = pend-pstart;
		//std::cout << "plen = " << plen << std::endl;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		//std::cout << "to_match = " << to_match << std::endl;
		set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
                                            		    pend-to_match, to_match);
		std::cout << "search = " << search << std::endl;

		uint_t r = ef.rank1(search+1);
		std::cout << "r = " << r << std::endl;
		uint_t selected;
		uint_t s1 = ef.select1(r-1);
		std::cout << "s1 = " << s1 << std::endl;

		//exit(1);

		//ATAA
		//ATAAC

		uint_t mlen = 0;
		if(s1 == search)
		{
			if(plen <= this->len)
			{
				//std::cout << " 1 return " << Suff[bv.select(r-1)] << " " << to_match << " false" << std::endl;
				//exit(1);
				return std::make_tuple(this->Suff[bv.select(r-1)],to_match,false);
			}
			else{ /*std::cout << "qui?" << std::endl;*/ selected = s1; r--; mlen = to_match; }
		}
		else if(r == ef.no_ones())
		{
			//std::cout << "quaaaa?" << std::endl;
			//std::cout << __builtin_clz(search ^ s1) << std::endl;
			//std::cout << ((sizeof(uint_t)*8)-(this->len*2)) << std::endl;
			//std::cout << (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2 << " " << s1 << " " << r-- << std::endl;

			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2; selected = s1; r--;
		}
		else
		{
			//std::cout << "quIAIAIA?" << std::endl;
			uint_t s2 = ef.select1(r);
			//std::cout << "s2 = " << s2 << std::endl;
			int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);
			//std::cout << m1 << " " << m2 << std::endl;
			if( m1 > m2 )
			{
				//std::cout << "S1 = " << s1 << std::endl;
				selected = s1;
				r--;
				mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
			}
			else
			{
				//std::cout << "S2 = " << s2 << std::endl;
				selected = s2;
				mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
			}
		}

		//std::cout << "SEARCH = " << search << std::endl;

		//std::cout << "mlen = " << mlen << std::endl;
		//std::cout << "selected = " << selected << std::endl;
		//std::cout << "r = " << r << std::endl;

		//if((mlen < to_match) or (plen == to_match)){
		if(mlen < to_match){
			//std::cout << "--->" << this->Suff[bv.select(r)] << " " << to_match << " true" << std::endl;
			//exit(1);
			return std::make_tuple(this->Suff[bv.select(r)],mlen,true);
		}
		else if(plen == to_match){ return std::make_tuple(this->Suff[bv.select(r)],to_match,false); }

		// else run binary search

		uint_t low, mid, high;
		int_t lcp_low, lcp_high, lcp_mid;
		low  = bv.select(r);
		high = bv.select(r+1);

		// std::cout << low << " - " << high << std::endl;

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
			////std::cout << "computing LCS with : " << this->Suff[mid] << std::endl;			
			auto j = O->LCS_char(pattern,pend-1,this->Suff[mid]); 
	
			if(j.first == plen)
				return std::make_tuple(this->Suff[mid],plen,false); 		

			////std::cout << "pend: " << pend << std::endl;
			////std::cout << "-- " << j.first << " - " << j.second << " <-> " << pattern[pend-j.first-1] << std::endl; 
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

		////std::cout << "lcp_low= " <<  lcp_low << " - lcp_high = " << lcp_high << std::endl;

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

		//std::cout << "mid: " << mid << " -> " << this->Suff[mid] << std::endl;
		//std::cout << "lcp_mid: " << lcp_mid << std::endl;

		//std::cout << " 2 return " << this->Suff[mid] << " " << lcp_mid << " " << (lcp_mid != plen) << std::endl;
		return std::make_tuple(this->Suff[mid],lcp_mid,(lcp_mid != plen));

		//exit(1);
		/*
		////std::cout << "locate_longest_prefix = " << pattern.substr(pstart,pend-pstart) << std::endl;
		// initialize binary search parameters
		uint_t low, mid, high, plen;
		int_t lcp_low, lcp_high, lcp_mid;
		plen = pend - pstart;
		//low  = this->alph[dna_to_code_table[pattern[pend-1]]];
		//high = this->alph[dna_to_code_table[pattern[pend-1]]+1];
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
		////std::cout << "low = " << low << " - " << " high = " << high << std::endl;

		while( high-low > 1 )
		{
			////std::cout << "computing LCS with : " << this->Suff[mid] << std::endl;			
			auto j = O->LCS_char(pattern,pend-1,this->Suff[mid]); 
	
			if(j.first == plen)
				return std::make_tuple(this->Suff[mid],plen,false); 		

			////std::cout << "pend: " << pend << std::endl;
			////std::cout << "-- " << j.first << " - " << j.second << " <-> " << pattern[pend-j.first-1] << std::endl; 
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

		////std::cout << "lcp_low= " <<  lcp_low << " - lcp_high = " << lcp_high << std::endl;

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
		*/
		// return std::make_tuple(0,1,false);
	}

private:
	//
	text_oracle* O;
	//
	usafe_t N;
	//
	int_t len;
	//
	sdsl::int_vector<> Suff;
	//
	elias_fano_ds ef;
	//
	bitvector bv;
};

}

#endif // SUFFIXIENT_ARRAY_BASELINE_HPP_