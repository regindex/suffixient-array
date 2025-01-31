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

void help(){

	cout << "suffixient [options]" << endl <<
	"Input: non-empty ASCII file without character 0x0, from standard input. Output: smallest suffixient-nexessary set." << endl <<
	"Warning: if 0x0 appears, the standard input is read only until the first occurrence of 0x0 (excluded)." << endl <<
	"Options:" << endl <<
	"-h          Print usage info." << endl <<
	"-l          Store LCS information between supermaximal substrings. Default: false." << endl << 
	"-o <arg>    Output files basepath. If not specified, output is streamed to standard output in human-readable format." << endl <<
	"-t          Remap original alphabet to a contiguos range 0,1,...,N. Default: false." << endl <<
	"-s          Sort output suffixient set. Default: false." << endl <<
	"-p          Print to standard output size of suffixient set. Default: false." << endl <<
	"-r          Print to standard output number of equal-letter runs in the BWT of reverse text. Default: false." << endl <<
	"-i          Print to standard output statistics on the LCS value distribution. Default: false." << endl;
	exit(0);
} 

inline void eval(uint64_t sigma, int64_t m, vector<lcp_maxima>& R,  
								 vector<uint64_t>& S, vector<int64_t>& L,
								 vector<uint64_t>& A, vector<int64_t>& last)
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

		    	if(last[c] != -1) L[last[c]] = max(L[last[c]],R[c].lcs + 1);
		    	last[c] = L.size()-1;
		    }
		    // update to inactive state
		    R[c] = {m,0,m,false};
		  }
}

int main(int argc, char** argv){

	if(argc < 2) help();

	string output_basepath;

	bool sort = false;
	bool rho = false;
	bool runs = false;
	bool remap = false;
	bool lcs_info = false;
	bool lcs_stats = false;

	int opt;
	while ((opt = getopt(argc, argv, "prshtlio:")) != -1){
		switch (opt){
			case 'h':
				help();
			break;
			case 'l':
				lcs_info=true;
			break;
			case 'o':
				output_basepath = string(optarg);
			break;
			case 's':
				sort=true;
			break;
			case 'p':
				rho=true;
			break;
			case 'r':
				runs=true;
			break;
			case 't':
				remap=true;
			break;
			case 'i':
				lcs_stats=true;
			break;
			default:
				help();
			return -1;
		}
	}
	
	cache_config cc;
	uint64_t N = 0; //including 0x0 terminator
	int sigma = 1; // alphabet size (including terminator 0x0)
	int64_t m = std::numeric_limits<int64_t>::max();
	uint64_t bwtruns = 1;

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

		if(not remap)
		{
			for(uint64_t i = 0; i < N - 1; ++i)
			{
				uint8_t c = in[N - i - 2];
				T[i] = c;
				sigma = max(sigma,int(c));
			}
			sigma++;
		}
		else
			for(uint64_t i = 0; i < N - 1; ++i)
			{
				uint8_t c = in[N - i - 2];
				if(char_to_int[c] == 0) char_to_int[c] = sigma++;
				T[i] = char_to_int[c];
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
	vector<uint64_t> S;
	vector<int64_t>  L;
	vector<uint64_t> A(sigma,0);
	vector<int64_t>  last(sigma,-1);

	for(uint64_t i=1;i<N;++i)
	{
		m = std::min(m,int64_t(LCP[i]));

		if(BWT(i) != BWT(i-1))
		{
			eval(sigma,m,R,S,L,A,last);

			for(uint64_t ip = i-1; ip < i+1; ++ip)
				if(int64_t(LCP[i]) > R[BWT(ip)].len)
					R[BWT(ip)] = {int64_t(LCP[i]),N - SA[ip],R[BWT(ip)].lcs,true}; 
      // reset LCP value
      m = std::numeric_limits<int64_t>::max();
      // increment number of runs
      bwtruns++;
		}
	}
  // evaluate last active candidates
  eval(sigma,-1,R,S,L,A,last); 

  // remove chached files
  sdsl::remove(cache_file_name(conf::KEY_TEXT,cc));
  sdsl::remove(cache_file_name(conf::KEY_SA,  cc));
  sdsl::remove(cache_file_name(conf::KEY_ISA, cc));
  sdsl::remove(cache_file_name(conf::KEY_LCP, cc));

  if(sort) std::sort(S.begin(),S.end());

  if(output_basepath.length()==0)
  {
  	for(auto x:S) cout << x << " ";
  	cout << endl;
  	if(lcs_info)
  	{
	  	for(auto x:L) cout << x << " ";
	  	cout << endl;
		}
  }
  else
  {
  	string suffixient_output_basepath = output_basepath + ".suff";
    ofstream ofs(suffixient_output_basepath, ios::binary);
    for (const auto& x : S) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
    ofs.close();
  	string suff_alph_output_basepath = output_basepath + ".alph";
  	ofs = ofstream(suff_alph_output_basepath, ios::binary);
  	for (const auto& x : A) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
  	ofs.close();
  	if(lcs_info)
  	{
	  	string lcs_output_basepath = output_basepath + ".lcs";
	  	ofs = ofstream(lcs_output_basepath, ios::binary);
	  	for (const auto& x : L) { ofs.write(reinterpret_cast<const char*>(&x), STORE_SIZE); }
	  	ofs.close();
		}
  }

  if(rho) cout << "Size of smallest suffixient set: " << S.size() << endl;
  if(runs) cout << "Number of equal-letter BWT(rev(T)) runs: " << bwtruns << endl;
  if(lcs_stats)
  {
  	size_t N = 34;
  	std::vector<uint64_t> freq(N,0);
  	for(size_t i=0;i<L.size();++i)
  	{
  		if(L[i] <= 32)
  			freq[L[i]]++;
  		else
  			freq[N-1]++;
  	}

  	for(size_t i=1;i<N;++i){ freq[i] += freq[i-1]; }
  	for(size_t i=1;i<N;++i){ std::cout << "i: " << i << " - " 
  		                                 << (double(freq[i])/L.size())*100 << "%" << std::endl; }
  }
	
  return 0;
}