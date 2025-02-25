// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <sdsl/construct.hpp>
#include <set>
#include <limits>
#include <algorithm>

using namespace std;
using namespace sdsl;

#define STORE_SIZE 5

int_vector<8> T;
int_vector_buffer<> SA;
int_vector_buffer<> LCP;

cache_config cc;
uint64_t N = 0; //including 0x0 terminator
int sigma = 1; // alphabet size (including terminator 0x0)
int64_t m = std::numeric_limits<int64_t>::max();

inline uint8_t BWT(uint64_t i){ return SA[i] == 0 ? 0 : T[SA[i] - 1]; }

struct lcp_maxima
{
	int64_t len;
	uint64_t pos;
	int64_t lcs;
	bool active;
	lcp_maxima(int64_t len_, uint64_t pos_, int64_t lcs_, bool active_) :
					   len(len_), pos(pos_), lcs(lcs_), active(active_) {}
};

void help()
{
	cout << "Software to compute all ingredients necessary to construct different variants for the suffixient index." << endl <<
	        "one-pass-index [options]" << endl <<
	"Input: non-empty ASCII file without character 0x0, from standard input." << endl <<
	"Output: depending on the index variant: smallest suffixient vector, supermaximal substrings LCS vector, right-extension multiplicities, Prefix array. " <<
			    "All integers are store in binary format using 5 bytes per entry." << endl <<
	"Warning: if 0x0 appears, the standard input is read only until the first occurrence of 0x0 (excluded)." << endl <<
	"Options:" << endl <<
	"-h          Print usage info." << endl <<
	"-t          Select the index variant for which you want to compute the components" <<
	             "(z-fastTrie|suffixientArray|prefixArray). Default: false." << endl << 
	"-o <arg>    Output files basepath. If not specified, output is streamed to standard output in human-readable format." << endl;
	exit(0);
} 

inline void eval(uint64_t sigma, int64_t m, vector<lcp_maxima>& R,  
			     	   vector<uint64_t>& S,    vector<int64_t>& Len, 
			     							   vector<int64_t>& last)
{
	for(uint8_t c = 1; c < sigma; ++c)
	  if(m < R[c].len)
		  {
		    // process an active candidate
		    if(R[c].active)
		    {
		    	S.push_back(R[c].pos - 1);
		    	Len.push_back(R[c].lcs + 1);

		    	if(last[c] != -1) Len[last[c]] = max(Len[last[c]],R[c].lcs + 1);
		    	last[c] = Len.size()-1;
		    }
		    // update to inactive state
		    R[c] = {m,0,m,false};
		  }
}

inline void eval(uint64_t sigma, int64_t m, vector<lcp_maxima>& R,  
				       vector<uint64_t>& S,    vector<int64_t>& L, 
				                              vector<uint64_t>& A)
{
	for(uint8_t c = 1; c < sigma; ++c)
	  if(m < R[c].len)
		  {
		    // process an active candidate
		    if(R[c].active)
		    {
		    	S.push_back(R[c].pos - 1);
		    	L.push_back(R[c].lcs + 1);
		    	A[c]++;
		    }
		    // update to inactive state
		    R[c] = {m,0,m,false};
		  }
}

inline void one_pass_algorithm(vector<lcp_maxima>& R, vector<uint64_t>& S, 
	                           vector<int64_t>&  Len, vector<int64_t>&  last)
{
	for(uint64_t i=1;i<N;++i)
	{
		m = std::min(m,int64_t(LCP[i]));

		if(BWT(i) != BWT(i-1))
		{
			eval(sigma,m,R,S,Len,last);

			for(uint64_t ip = i-1; ip < i+1; ++ip)
				if(int64_t(LCP[i]) > R[BWT(ip)].len)
					R[BWT(ip)] = {int64_t(LCP[i]),N - SA[ip],R[BWT(ip)].lcs,true}; 

			// reset LCP value
			m = std::numeric_limits<int64_t>::max();
		}
	}
	// evaluate last active candidates
	eval(sigma,-1,R,S,Len,last); 
}

inline void one_pass_algorithm(vector<lcp_maxima>& R, vector<uint64_t>& S, 
	                           vector<int64_t>&    L, vector<uint64_t>& A)
{
	for(uint64_t i=1;i<N;++i)
	{
		m = std::min(m,int64_t(LCP[i]));

		if(BWT(i) != BWT(i-1))
		{
			eval(sigma,m,R,S,L,A);

			for(uint64_t ip = i-1; ip < i+1; ++ip)
				if(int64_t(LCP[i]) > R[BWT(ip)].len)
					R[BWT(ip)] = {int64_t(LCP[i]),N - SA[ip],R[BWT(ip)].lcs,true}; 

			// reset LCP value
			m = std::numeric_limits<int64_t>::max();
		}
	}
	// evaluate last active candidates
	eval(sigma,-1,R,S,L,A); 
}

int main(int argc, char** argv)
{
	if(argc < 2) help();

	string output_basepath, index_variant;

	int opt;
	while ((opt = getopt(argc, argv, "ht:o:")) != -1)
	{
		switch (opt)
		{
			case 'h':
				help();
			break;
			case 't':
				index_variant = string(optarg);
			break;
			case 'o':
				output_basepath = string(optarg);
			break;
			default:
				help();
			return -1;
		}
	}

	if(index_variant != "z-fastTrie" and index_variant != "suffixientArray"
		                             and index_variant !=     "prefixArray")
	{
		std::cerr << "Error! select a valid index variant..." << std::endl;
		help();
	}

	{
		string in;
		getline(cin,in,char(0));
		N = in.size() + 1;

		if(N<2){
			cerr << "Error: empty text" <<  endl;
			help();
		}

		T = int_vector<8>(N - 1);

		vector<uint8_t> char_to_int(256, 0); //map chars to 0...sigma-1. 0 is reserved for term.

		{
			for(uint64_t i = 0; i < N - 1; ++i)
			{
				uint8_t c = in[N - i - 2];
				T[i] = c;
				sigma = max(sigma,int(c));
			}
			sigma++;
		}

		append_zero_symbol(T);
		store_to_cache(T, conf::KEY_TEXT, cc);
		construct_sa<8>(cc);
		construct_lcp_kasai<8>(cc);
		SA = int_vector_buffer<>(cache_file_name(conf::KEY_SA, cc));
		LCP = int_vector_buffer<>(cache_file_name(conf::KEY_LCP, cc));
	}

	//vector with candidate suffixient right-extensions
	vector<lcp_maxima> R(sigma,{-1,0,-1,false}); 
	vector<uint64_t> S, A(sigma,0);
	vector<int64_t>  L, Len, last(sigma,-1);

	if(index_variant == "z-fastTrie")
	{
		one_pass_algorithm(R,S,Len,last);
		{
			string suffixient_output_basepath = output_basepath + ".suff";
		  	ofstream ofs(suffixient_output_basepath, ios::binary);
		  	for (const auto& x : S) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
		  	ofs.close();
		}
		{
		  	string lens_output_basepath = output_basepath + ".lens";
		  	ofstream ofs = ofstream(lens_output_basepath, ios::binary);
		  	for (const auto& x : Len) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
		  	ofs.close();
		}
	}
	else if(index_variant == "suffixientArray")
	{
		one_pass_algorithm(R,S,L,A);
		{
			string suffixient_output_basepath = output_basepath + ".suff";
		  	ofstream ofs(suffixient_output_basepath, ios::binary);
		  	for (const auto& x : S) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
		  	ofs.close();
		}
		{
		  	string lcs_output_basepath = output_basepath + ".lcs";
		  	ofstream ofs = ofstream(lcs_output_basepath, ios::binary);
		  	for (const auto& x : L) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
		  	ofs.close();
		}
		{
	  	string suff_alph_output_basepath = output_basepath + ".mult";
	  	ofstream ofs = ofstream(suff_alph_output_basepath, ios::binary);
	  	for (const auto& x : A) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
	  	ofs.close();
		}
	}
	else if(index_variant == "prefixArray")
	{
		std::cout << "N = " << N << std::endl;
		std::cout << "STORE SIZE = " << STORE_SIZE << std::endl;
	  	string pa_output_basepath = output_basepath + ".pa";
	  	ofstream ofs = ofstream(pa_output_basepath, ios::binary);
	  	for (size_t i=1;i<N;++i)
	  		{ size_t x = N - SA[i] - 1; ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
	  	ofs.close();
	}

	// remove chached files
	sdsl::remove(cache_file_name(conf::KEY_TEXT,cc));
	sdsl::remove(cache_file_name(conf::KEY_SA,  cc));
	sdsl::remove(cache_file_name(conf::KEY_ISA, cc));
	sdsl::remove(cache_file_name(conf::KEY_LCP, cc));
	
	return 0;
}