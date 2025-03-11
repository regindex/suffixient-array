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
		                             int_t len_, bool_t verbose = true)
	{
		usafe_t a = 0, i = 0, l = 0; 
	    std::vector<uchar_t> buffer(STORE_SIZE,0);
		{
			O = O_;
		}
		this->len = len_;
		if(this->len > 16 and sizeof(uint_t) == 4)
		{ 
			std::cerr << "The stored substrings length must be smaller or equal than 16 when using words of 4 bytes!" << std::endl; 
			std::cerr << "Please, run the 64 bits executable..." << std::endl; 
			exit(1);
		}
		if(this->len > 32)
		{ 
			std::cerr << "The supermaximal substrings length must be smaller or equal than 32!" << std::endl; 
			std::cerr << "Please, use -l flag to select a smaller additional space percentage..." << std::endl; 
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

					  	set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&offset),text_buffer);
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

			set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&offset),text_buffer);
			onset.push_back(offset);

			file_text.close();

			this->bv = bitvector(tmp_bv);
			this->ef = elias_fano_ds(onset,pow(SIGMA_DNA,this->len));

			{
				usafe_t bv_density = 0;
				for(usafe_t i=0;i<tmp_bv.size();++i){ if(tmp_bv[i]){ bv_density++; } }
				std::cout << "bv_density: " << (double)bv_density/tmp_bv.size() << std::endl;
				std::cout << "ef_density: " << (double)onset.size()/pow(SIGMA_DNA,this->len) << std::endl;
			}

			file_text.close();
		}

	    file_suff.close();
	    file_lcs.close();
	}

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
	/*
	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix(std::string& pattern,uint_t pstart,uint_t pend) const
	{
		uint_t search = 0;
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint_t mlen = 0;

		std::cout << "Heuristic len: " << this->len << std::endl;
		set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
                                            		    pend-to_match, to_match);
		//std::cout << "search = " << search << " - " << pstart << " - " << pend << std::endl;
		uint_t r = ef.rank1(search+1);

		//std::cout << "r " << r << std::endl;

		if(r == 0)
		{
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else if(r == ef.no_ones())
		{
			r--;
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else
		{
			uint_t s1 = ef.select1(r-1);
			if(search == s1){ r--; mlen = to_match; }
			else
			{
				uint_t s2 = ef.select1(r);
				int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);

				if(m1 > m2)
				{
					r--;
					mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
				}
				else{ mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2; }
			}
		}
		
		//else
		//{
		//	uint_t s1 = ef.select1(r-1), s2 = ef.select1(r);
		//	int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);
		//
		//	if(m1 > m2)
		//	{
		//		r--;
		//		mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
		//	}
		//	else{ mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2; }
		//}

		if(mlen < to_match){ return std::make_tuple(this->Suff[bv.select(r)],mlen,true); }
		else if(plen == to_match){ return std::make_tuple(this->Suff[bv.select(r)],to_match,false); }
		else
		{
			uint_t low, mid, high;
			int_t lcp_low, lcp_high, lcp_mid;
			low  = bv.select(r);
			high = bv.select(r+1);

			//std::cout << "low " << low << " high " << high << std::endl;

			// stop if first pattern character doesn't occur in the text
			if((high - low) > 0)
			{ 
				// --> devo lavorare qua dentro!!!
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
	}
	*/

	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix(std::string& pattern,uint_t pstart,uint_t pend) const
	{
		uint_t search = 0;
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint_t mlen = 0;

		//std::cout << "Heuristic len: " << this->len << std::endl;
		set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
                                            		    pend-to_match, to_match);
		//std::cout << "search = " << search << " - " << pstart << " - " << pend << std::endl;
		uint_t r = ef.rank1(search+1);

		//std::cout << "r " << r << std::endl;

		if(r == 0)
		{
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else if(r == ef.no_ones())
		{
			r--;
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else
		{
			uint_t s1 = ef.select1(r-1);
			if(search == s1){ r--; mlen = to_match; }
			else
			{
				uint_t s2 = ef.select1(r);
				int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);

				if(m1 > m2)
				{
					r--;
					mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
				}
				else{ mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2; }
			}
		}

		if(mlen < to_match){ return std::make_tuple(this->Suff[bv.select(r)],mlen,true); }
		else if(plen == to_match){ return std::make_tuple(this->Suff[bv.select(r)],to_match,false); }
		else
		{
			uint_t low, mid, high;
			int_t lcp_low, lcp_high, lcp_mid;
			low  = bv.select(r);
			high = bv.select(r+1);

			//std::cout << "low " << low << " high " << high << std::endl;

			// stop if first pattern character doesn't occur in the text
			if((high - low) > 0)
			{ 
				// --> devo lavorare qua dentro!!!
				high--;
				lcp_low = lcp_high = -1; 
				mid = (low+high)/2;
				if(plen == 1)
					return std::make_tuple(this->Suff[mid],1,false);
			}
			else
				return std::make_tuple(-1,0,true);


			if( high-low == 0 )
			{
				lcp_mid  = O->LCS_char(pattern,pend-1,Suff[mid]).first;
				return std::make_tuple(this->Suff[mid],lcp_mid,(lcp_mid != plen));
			}
			else
			{
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
		}
	}

	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix_opt(std::string& pattern,uint_t pstart,uint_t pend) const
	{
		/*
		std::cout << "PROVETTA" << std::endl;
		std::cout << "no ones ef = " << ef.no_ones() << std::endl;
		usafe_t len_1_int = 0;
		usafe_t len_tot = 0;
		for(usafe_t i=0;i<ef.no_ones();++i)
		{
			usafe_t low  = bv.select(i), high = bv.select(i+1);
			if((high - low) == 1){ len_1_int++; }
			len_tot += (high - low);
		}
		std::cout << "no len_1_int = " << len_1_int << std::endl;
		std::cout << "len tot = " << len_tot << std::endl;
		usafe_t no_ones = 0;
		for(usafe_t i=0;i<bv.size();++i)
		{
			if(static_cast<bool>(bv[i])){ no_ones++; }
		}
		std::cout << "size = " << bv.size() << " no_ones = " << no_ones << std::endl;

		exit(1);
		*/
		uint_t search = 0;
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint_t mlen = 0;

		set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
                                            		    pend-to_match, to_match);
		uint_t r = ef.rank1(search+1);

		if(r == 0)
		{
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else if(r == ef.no_ones())
		{
			r--;
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else
		{
			uint_t s1 = ef.select1(r-1);
			if(search == s1){ r--; mlen = to_match; }
			else
			{
				uint_t s2 = ef.select1(r);
				int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);

				if(m1 > m2)
				{
					r--;
					mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
				}
				else{ mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2; }
			}
		}

		if(mlen < to_match)      { return std::make_tuple(this->Suff[bv.select(r)],mlen,true);      }
		else if(plen == to_match){ return std::make_tuple(this->Suff[bv.select(r)],to_match,false); }
		else
		{
			uint_t low, mid, high;
			int_t lcp_low, lcp_high, lcp_mid;
			low  = bv.select(r); high = bv.select(r+1);

			// stop if first pattern character doesn't occur in the text
			if((high - low) > 0)
			{ 
				high--;
				lcp_low = lcp_high  = -1; 
				mid = (low+high)/2; 
				if(plen == 1){ return std::make_tuple(this->Suff[mid],1,false); }
			}
			else{ return std::make_tuple(-1,0,true); }

			if((high - low) == 0) { return std::make_tuple(this->Suff[mid],mlen,false); }
			else{
				while( high-low > 1 )
				{			
					auto j = O->LCS_char(pattern,pend-1,this->Suff[mid]); 
			
					if(j.first == plen){ return std::make_tuple(this->Suff[mid],plen,false); } 		

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
		}
	}
	/*
	std::tuple<uint_t,uint_t,bool_t> 
		locate_longest_prefix(sdsl::int_vector<2>& pattern,uint_t pstart,uint_t pend) const
	{
		uint_t search = 0;
		uint_t plen = pend-pstart;
		uint_t to_match = std::min(static_cast<uint_t>(this->len),plen);
		uint_t mlen = 0;


		//set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
        //                                    		    pend-to_match, to_match);

		search = pattern.get_int((pend-to_match)*2,to_match*2);
		//if( ((this->len - to_match)*2) > 0  )
			search <<= (this->len - to_match)*2;

		//std::cout << "search = " << search << " - " << pstart << " - " << pend << std::endl;
		//exit(1);

		//set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&search), pattern, this->len, 
        //                                    		    pend-to_match, to_match);
		uint_t r = ef.rank1(search+1);

		//std::cout << "r " << r << std::endl;

		if(r == 0)
		{
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else if(r == ef.no_ones())
		{
			r--;
			uint_t s1 = ef.select1(r);
			mlen = (__builtin_clz(search ^ s1)-((sizeof(uint_t)*8)-(this->len*2)))/2;
		}
		else
		{
			uint_t s1 = ef.select1(r-1);
			if(search == s1){ r--; mlen = to_match; }
			else
			{
				uint_t s2 = ef.select1(r);
				int_t m1 = __builtin_clz(search ^ s1), m2 = __builtin_clz(search ^ s2);

				if(m1 > m2)
				{
					r--;
					mlen = (m1 - ((sizeof(uint_t)*8)-(this->len*2)))/2;
				}
				else{ mlen = (m2 - ((sizeof(uint_t)*8)-(this->len*2)))/2; }
			}
		}

		if(mlen < to_match){ return std::make_tuple(this->Suff[bv.select(r)],mlen,true); }
		else if(plen == to_match){ return std::make_tuple(this->Suff[bv.select(r)],to_match,false); }
		else
		{
			//std::cout << "ENTRA NELLA BINARY SEARCH" << std::endl;
			uint_t low, mid, high;
			int_t lcp_low, lcp_high, lcp_mid;
			low  = bv.select(r);
			high = bv.select(r+1);

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

			//std::cout << "low " << low << " high " << high << std::endl;

			while( high-low > 1 )
			{			
				auto j = O->LCS_char(pattern,pend-1,this->Suff[mid]);
				//std::cout << "--> " << j.first << " | " << j.second << std::endl; 
		
				if(j.first == plen)
					return std::make_tuple(this->Suff[mid],plen,false); 		

				if(j.second > code_to_dna_table[pattern[pend-j.first-1]]){
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
	}
	*/
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