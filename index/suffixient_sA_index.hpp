// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_index: DESCRIPTION
 */

#ifndef SUFFIXIENT_sA_INDEX_HPP_
#define SUFFIXIENT_sA_INDEX_HPP_

#include <chrono>
#include <malloc_count.h>
// suffixient array data structure
#include <suffixient_array_baseline.hpp>
// suffixient array elias-fano augmented data structure 
#include <suffixient_array_elias_fano.hpp>
// suffix array data structure
#include <suffix_array_binary_search.hpp>
// Text oracle data structure
//#include <LZ77_compressed_text.hpp>
#include <LZ77_compressed_text-modified.hpp>
//#include <static_compressed_text.h>
// K-order entropy compressed text oracle
#include <Hk_Ferragina_Venturini.hpp>
// Uncompressed text oracle data structure
#include <uncompressed_text_oracle.hpp>
// Bitpacked text oracle data structure
#include <bitpacked_text_oracle.hpp>
// Relative Lempel Ziv compressed text
//#include <RLZ_DNA.hpp>
#include <RLZ_DNA_o2.hpp>

namespace suffixient{

template<class suffArray, class textOracle>
class suffixient_sA_index{

public:
	// empty constructor
	suffixient_sA_index(){}

	// build suffixient index by indexing the supermaximal extensions
	void build(const std::string &text_filepath, const uint8_t stored_len,
		 									     const bool_t verbose = false)
	{
		std::cout << "Step 1) Constructing text oracle data structure for " 
		          << text_filepath << std::endl;
	  	O.build(text_filepath);
		std::cout << "Step 2) Constructing the prefix search data structure for " 
				  << text_filepath+".suff" << std::endl;
		S.build(text_filepath, &O, stored_len, verbose);
	}

	usafe_t store(const std::string &index_filepath)
	{
		std::ofstream out(index_filepath);

		usafe_t w_bytes = 0;
		w_bytes += O.size();
		w_bytes += S.store(out);
		
		out.close();

		return w_bytes;
	}

	void load(const std::string &oracle_filepath, 
		      const std::string &index_filepath)
	{
		O.load(oracle_filepath);

		std::ifstream in(index_filepath);
		S.load(in,&O);
		in.close();
	}

	void find_MEMs(string& pattern, std::ofstream& output)
	{
		////std::cout << "############Finding MEMs for = " << pattern << std::endl;
		size_t i = 0, l = 0, m = pattern.size();
		int64_t pstart = 0;
		while(i < m)
		{
			////std::cout << "i = " << i << std::endl;
			////std::cout << "pstart: " << pstart << " pend: " << (i+1) << std::endl;
			////std::cout << "right maximal substring = ";
			/*for(int k=pstart;k<(i+1);++k)
				std::cout << pattern[k];
			std::cout << std::endl;*/

			//auto j = this->S.locate_longest_prefix(pattern,pstart,(i+1)-pstart);
			////std::cout << "locating prefix: " << pstart << " " << i+1 << std::endl;
			auto j = this->S.locate_longest_prefix(pattern,pstart,i+1);
			////std::cout << "finish locate" << std::endl;
			////std::cout << std::get<0>(j) << " " << std::get<1>(j) << " " << std::get<2>(j) << std::endl;
			size_t b = std::get<1>(j);
			// if the current right-extension does not occur in the text
			if(b == 0)
			{
				if(l > 0)
				{
					std::cout << "(" << i-l << "," << l << ") ";
					output << "(" << i-l << "," << l << ") ";
				}
				std::cout << "(" << i << "," << b << ") ";
				output << "(" << i << "," << b << ") ";
				pstart = i = i+1;
				l = 0;
			}
			else
			{
				if(b <= l)
				{
					std::cout << "(" << i-l << "," << l << ") ";
					output << "(" << i-l << "," << l << ") ";
					pstart = i-l+1;
					////std::cout << "pstart: " << pstart << std::endl;
				}

				size_t f = O.LCP(pattern,i+1,std::get<0>(j)+1);
				////std::cout << "f: " << f << std::endl;
				i = i + f + 1;
				////std::cout << "i: " << i << std::endl;
				l = b + f;
				////std::cout << "l: " << l << std::endl;
			}
		}
		if(l > 0)
		{
		std::cout << "(" << i-l << "," << l << ")" << std::endl;
		output << "(" << i-l << "," << l << ")" << std::endl;
		}
		else
		{
			std::cout << std::endl;
			output << std::endl;
		}
	}
	/*
	template<typename stype>
	std::pair<usafe_t,double> locate_one_occurrence(stype pattern)
	{
		//std::cout << "#### LOCATING " << pattern << std::endl;
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, occ = -1, m = pattern.size();
		bool_t mismatch_found;
		while(i-1 < m)
		{
			//std::cout << "i: " << i << std::endl;
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			//std::cout << "--i: " << i << std::endl;
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return std::make_pair(-1,0);

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(occ-m+1,duration.count());
	}
	template<typename stype>
	std::pair<safe_t,double> locate_one_occurrence_heuristic(stype pattern)
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t occ = -1, m = pattern.size(), 
		                  i = std::min(S.get_len(),static_cast<int_t>(m));
		bool_t mismatch_found;

		while(i > 0)
		{
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			mismatch_found = std::get<2>(j);
			if(not mismatch_found)
			{
				occ = std::get<0>(j);
				if(i < m)
				{
					usafe_t f = O.LCP(pattern,i,occ+1);
					i = i + f + 1;
					occ = occ + f;
				}
				break;
			}
			i--;
		}

		while(i-1 < m)
		{
			//std::cout << "I= " << i << std::endl;
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			//std::cout << std::get<0>(j) << " " << std::get<1>(j) << " " << std::get<2>(j) << std::endl;
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return std::make_pair(-1,0);

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(occ-m+1,duration.count());
	}
	template<typename stype>
	std::pair<usafe_t,double> locate_longest_prefix(stype prefix)
	{
		auto start = std::chrono::high_resolution_clock::now();

		auto j = this->S.locate_longest_prefix(prefix,0,prefix.size());

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(std::get<0>(j)-prefix.size()+1,duration.count());
	}
	*/

	std::pair<usafe_t,double> locate_one_occurrence(std::string pattern)
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t i = 1, occ = -1, m = pattern.size();
		bool_t mismatch_found;
		while(i-1 < m)
		{
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return std::make_pair(-1,0);

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(occ-m+1,duration.count());
	}

	std::pair<safe_t,double> locate_one_occurrence_heuristic(std::string pattern)
	{
		auto start = std::chrono::high_resolution_clock::now();

		usafe_t occ = -1, m = pattern.size(), 
		                  i = std::min(S.get_len(),static_cast<int_t>(m));
		bool_t mismatch_found, final_check = false;

		while(i > 0)
		{
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			mismatch_found = std::get<2>(j);
			if(not mismatch_found)
			{
				occ = std::get<0>(j);
				if(i < m)
				{
					usafe_t f = O.LCP(pattern,i,occ+1);
					i = i + f + 1;
					occ = occ + f;
				}
				break;
			}
			i--;
		}

		while(i-1 < m)
		{
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			mismatch_found = std::get<2>(j);

			if(mismatch_found)
				return std::make_pair(-1,0);

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
		}

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(occ-m+1,duration.count());
	}

	std::pair<safe_t,double> locate_one_occurrence_heuristic_ef_opt(std::string pattern)
	{
		auto start = std::chrono::high_resolution_clock::now();
		//std::cout << "Pattern = " << pattern << std::endl;

		usafe_t occ = -1, m = pattern.size(), 
		                  i = std::min(S.get_len(),static_cast<int_t>(m));
		bool_t mismatch_found;
		usafe_t matched_len;
		safe_t last_long_prefix = -1, last_occ = -1;
		//std::cout << "i = " << int(i) << std::endl;

		while(i > 0)
		{
			auto j = this->S.locate_longest_prefix(pattern,0,i);
			mismatch_found = std::get<2>(j);
			//std::cout << "mismatch_found = " << int(mismatch_found) << " - " << std::get<0>(j) << std::endl;
			if(not mismatch_found)
			{
				occ = std::get<0>(j);
				if(i < m)
				{
					usafe_t f = O.LCP(pattern,i,occ+1);
					i = i + f + 1;
					occ = occ + f;
				}
				else
				{
					std::chrono::duration<double> duration = 
							std::chrono::high_resolution_clock::now() - start;
					return std::make_pair(occ-m+1,duration.count());
				}
				break;
			}
			i--;
		}

		//std::cout << "i = " << int(i) << std::endl;

		while(i-1 < m)
		{
			auto j = this->S.locate_longest_prefix_opt(pattern,0,i);
			matched_len = std::get<1>(j); mismatch_found = std::get<2>(j);

			//std::cout << i << " - " << std::get<0>(j) << " : " << matched_len << " - " << mismatch_found << std::endl;

			if(mismatch_found){ return std::make_pair(-1,0); }
			if(matched_len < i){ last_long_prefix = i; last_occ = std::get<0>(j); }

			occ = std::get<0>(j);
			usafe_t f = O.LCP(pattern,i,occ+1);
			i = i + f + 1;
			occ = occ + f;
			//std::cout << "extended i = " << i << "extended occ = " << occ << std::endl;
		}

		if(last_long_prefix > -1)
		{
			//std::cout << "# " << last_long_prefix << " - " << last_occ << std::endl;
			usafe_t f = O.LCS(pattern,last_long_prefix-1,last_occ);
			//std::cout << "f = " << f << std::endl;
			if(f < last_long_prefix){ occ = -1; }
		}

		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		//exit(1);
		//std::cout << "return " << occ-m+1 << std::endl;

		return std::make_pair(occ-m+1,duration.count());
	}

	std::pair<usafe_t,double> locate_longest_prefix(std::string prefix)
	{
		auto start = std::chrono::high_resolution_clock::now();

		auto j = this->S.locate_longest_prefix(prefix,0,prefix.size());
		if( std::get<1>(j) != prefix.size() ) return std::make_pair(-1,0);
		
		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;

		return std::make_pair(std::get<0>(j)-prefix.size()+1,duration.count());
	}

	void locate_MEMs_fasta(std::string patternFile)
	{
		std::ifstream patterns(patternFile);
		std::ofstream output(patternFile+".mems");

		std::string line;
		size_t i=0,c=0;

		malloc_count_reset_peak();
		auto start = std::chrono::high_resolution_clock::now();

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				find_MEMs(line,output);
				c += line.size();
			}
			else{ output << line << std::endl; std::cout << line << std::endl; }
			i++;
		}

		patterns.close();
		output.close();

		output = std::ofstream(patternFile+".stats");
    
		std::chrono::duration<double> duration = 
				std::chrono::high_resolution_clock::now() - start;
		output << "Memory peak while searching for MEMs = " <<
				     malloc_count_peak() << " Kb" << std::endl;
		output << "Elapsed time while searching for MEMs = " <<
				     duration.count() << " sec" << std::endl;
		output << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl;
		output << "Elapsed time per pattern = " <<
				     (duration.count()/(i/2))*1000 << " milliSec" << std::endl;
		output << "Elapsed time per character = " <<
				     (duration.count()/(c))*1000000 << " microSec" << std::endl;

		output.close();

		return;
	}

	void run_exact_pattern_matching_fasta(
				 std::string patternFile, bool_t prefixArray  = false,
				                          bool_t runHeuristic = false)
	{
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".exactPM");

		std::string line, header;
		usafe_t i=0, c=0;
		std::pair<safe_t,double> o;
		double tot_duration = 0;

		malloc_count_reset_peak();

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				if(prefixArray){ o = locate_longest_prefix(line); }
				else if(runHeuristic)
					{ o = locate_one_occurrence_heuristic(line); }
				else{ o = locate_one_occurrence(line); }

				output << header << std::endl;
				if(o.first >= 0){ output << o.first << " " << line.size() << std::endl; }
				else{ output << "-1 " << line.size() << std::endl; }

				tot_duration += o.second;
				c += line.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000 << " milliSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000 << " microSec" << std::endl;
	}

	void run_exact_pattern_matching_fasta_ef_opt(std::string patternFile)
	{
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".exactPM");

		std::string line, header;
		usafe_t i=0, c=0;
		std::pair<safe_t,double> o;
		double tot_duration = 0;

		malloc_count_reset_peak();

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				o = locate_one_occurrence_heuristic_ef_opt(line);

				output << header << std::endl;
				if(o.first >= 0){ output << o.first << " " << line.size() << std::endl; }
				else{ output << "-1 " << line.size() << std::endl; }

				tot_duration += o.second;
				c += line.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000 << " milliSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000 << " microSec" << std::endl;
	}
	/*
	void run_exact_pattern_matching_fasta_bitpacked(
				 		  std::string patternFile, bool_t prefixArray  = false,
				                                   bool_t runHeuristic = false)
	{
		std::ifstream patterns(patternFile);
		std::ofstream   output(patternFile+".exactPM");

		std::string line, header;
		usafe_t i=0, c=0;
		std::pair<safe_t,double> o;
		double tot_duration = 0;

		malloc_count_reset_peak();

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				sdsl::int_vector<2> patt(line.size(),0);
				auto start = std::chrono::high_resolution_clock::now();

				for(usafe_t j=0;j<line.size();++j)
					{ patt[j] = dna_to_code_table[line[j]]; }

				std::chrono::duration<double> duration = 
					 std::chrono::high_resolution_clock::now() - start;
				tot_duration += duration.count();

				if(prefixArray){ o = locate_longest_prefix<sdsl::int_vector<2>>(patt); }
				else if(runHeuristic)
					{ o = locate_one_occurrence_heuristic<sdsl::int_vector<2>>(patt); }
				else{ o = locate_one_occurrence<sdsl::int_vector<2>>(patt); }
				//std::cout << "o = " << o.first << std::endl;
				output << header << std::endl;
				if(o.first >= 0){ output << o.first << " " << line.size() << std::endl; }
				else{ output << "-1 " << line.size() << std::endl; }

				tot_duration += o.second;
				c += line.size();
			}
			else{ header = line; }
			i++;
		}

		patterns.close();
		output.close();

		std::cout << "Memory peak while running pattern matching queries = " <<
				     malloc_count_peak() << " bytes" << std::endl
		          << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl
		          << "Number of patterns = " << i/2 
		 		  << ", Total number of characters = " << c << std::endl
		          << "Elapsed time per pattern = " <<
				     (tot_duration/(i/2))*1000 << " milliSec" << std::endl
		          << "Elapsed time per character = " <<
				     (tot_duration/(c))*1000000 << " microSec" << std::endl;
	}
	*/
	bool_t check_exact_pattern_matching_correctness(     std::string textFile,
						 std::string patternFile, bool_t prefixArray  = false,
						 						  bool_t runHeuristic = false)
	{
		std::cout << "Check correctness of exact pattern matching algorithm" << std::endl;
		if(prefixArray)
			std::cout << "Running binary search on the prefix array..." << std::endl;
		else
			std::cout << "Running binary search on the suffixient array..." << std::endl;
		if(runHeuristic)
			std::cout << "Running first prefix search heuristic..." << std::endl;

		suffixient::uncompressed_text_oracle T;
		std::cout << "Loading text in: " << textFile << std::endl;
		T.build(textFile);

		std::ifstream patterns(patternFile);

		std::string line;
		safe_t i=0;
		std::pair<safe_t,double> o;
		double tot_duration = 0;

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				//std::cout << line << " ";
				if(prefixArray){ o = locate_longest_prefix(line); }
				else if(runHeuristic)
					{ o = locate_one_occurrence_heuristic(line); }
				else{ o = locate_one_occurrence(line); }
				//std::cout << o << std::endl;
				tot_duration += o.second;

				if(o.first <= -1)
				{
					std::cerr << "Pattern: " << line 
					          << " not found!" << std::endl;
					return false;
				}
				else
				{
					std::string text_sub = T.display(o.first,o.first+line.size()-1);
					//std::cout << text_sub << std::endl;
					//exit(1);
					for(usafe_t j=0;j<line.size();++j)
						if(text_sub[j] != line[j])
						{
							std::cerr << "Error in the located pattern: "
							          << line << std::endl;
							return false;
						}
				}
			}
			i++;   
		}

		std::cout << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl;

		patterns.close();

		std::cout << "Everything's fine!" << std::endl;
		return true;
	}

	bool_t check_exact_pattern_matching_correctness_ef_opt(std::string textFile, std::string patternFile)
	{
		std::cout << "Check correctness of exact pattern matching algorithm" << std::endl;
		std::cout << "Running first prefix search heuristic with elias-fano optimization..." << std::endl;

		suffixient::uncompressed_text_oracle T;
		std::cout << "Loading text in: " << textFile << std::endl;
		T.build(textFile);

		std::ifstream patterns(patternFile);

		std::string line;
		safe_t i=0;
		std::pair<safe_t,double> o;
		double tot_duration = 0;

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				o = locate_one_occurrence_heuristic_ef_opt(line);
				tot_duration += o.second;

				if(o.first <= -1)
				{
					std::cerr << "Pattern: " << line 
					          << " not found!" << std::endl;
					return false;
				}
				else
				{
					std::string text_sub = T.display(o.first,o.first+line.size()-1);

					for(usafe_t j=0;j<line.size();++j)
						if(text_sub[j] != line[j])
						{
							std::cerr << "Error in the located pattern: "
							          << line << std::endl;
							return false;
						}
				}
			}
			i++;   
		}

		std::cout << "Elapsed time while running pattern matching queries = " <<
				     tot_duration << " sec" << std::endl;

		patterns.close();

		std::cout << "Everything's fine!" << std::endl;
		return true;
	}
	/*
	bool_t check_exact_pattern_matching_correctness_bitpacked(
						 std::string patternFile, bool_t prefixArray  = false,
						 						  bool_t runHeuristic = false)
	{
		std::ifstream patterns(patternFile);

		std::string line;
		safe_t i=0, o=0;

		while(std::getline(patterns, line))
		{
			if(i%2 != 0)
			{
				sdsl::int_vector<2> patt(line.size(),0);
				for(usafe_t j=0;j<line.size();++j)
					{ patt[j] = dna_to_code_table[line[j]]; }

				if(prefixArray)
					{ o = locate_longest_prefix<sdsl::int_vector<2>>(patt).first; }
				else if(runHeuristic)
					{ o = locate_one_occurrence_heuristic<sdsl::int_vector<2>>(patt).first; }
				else{ o = locate_one_occurrence<sdsl::int_vector<2>>(patt).first; }
				//std::cout << "o = " << o << std::endl;
				// check occurrence correctness
				if((o <= -1) or (O.LCP(patt,0,o) != patt.size()))
				{
					std::cerr << "Pattern: " << line 
					          << " not found!" << std::endl;
					return false;
				}
				//exit(1);
			}
			i++;
		}

		patterns.close();

		std::cout << "Everything's fine!" << std::endl;
		return true;
	}
	*/
private:
	// text oracle data structure
	textOracle O;
	// suffixient array search data structure
	suffArray S;
};

}

#endif
