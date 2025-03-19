// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_index: suffixient-index implementation based on z-fast trie
 */


#ifndef SUFFIXIENT_INDEX_HPP_
#define SUFFIXIENT_INDEX_HPP_

#include <chrono>
#include <malloc_count.h>
// common definitions
#include <common.hpp>
// prefix-search data structure
#include <CTriePP.hpp>
// LCP/LCS data structure
#include <LZ77_Kreft_Navarro_index.hpp>

namespace suffixient{

template<class prefixSearch, class lcpDS>
class suffixient_index{

public:
	// empty constructor
	suffixient_index(){}

	// build suffixient index by indexing the whole text prefixes
	std::pair<safe_t,safe_t> build(const std::string &text, const std::string &suffixient,
		                                                  bool verbose = false)
	{
		std::vector<uint_t> suffixient_set;
		safe_t tot_inserted_char = 0, tot_inserted_keywords = 0;
	    // Open the file in binary mode
	    std::ifstream file(suffixient, std::ios::binary);
	    if (!file.is_open()) {
	        std::cerr << "Error: Could not open file " << suffixient << std::endl;
	        exit(1);
	    }

	    safe_t number;
	    size_t count = 0;

	    while (file.read(reinterpret_cast<char*>(&number), 5))
	    {
	        assert(number > 0);
	        suffixient_set.push_back(number-1);
	    }
	    if (not file.eof()) { std::cerr << "An error occurred during reading." << std::endl; }

	    file.close();

	    std::sort(suffixient_set.begin(), suffixient_set.end());

	    file = std::ifstream(text, std::ios::binary);
	    if (!file.is_open()) {
	        std::cerr << "Error: Could not open file " << text << std::endl;
	        exit(1);
	    }

	    std::string last_prefix;
	    safe_t last_index = 0;
	    for(safe_t i=0;i<suffixient_set.size();++i)
	    {
	        std::string buffer(suffixient_set[i] - last_index + 1, '0');
	        file.read(&buffer[0], suffixient_set[i] - last_index + 1);
	        last_index = suffixient_set[i]+1;
	        
	        // invert the buffer
	        std::reverse(buffer.begin(),buffer.end());
	        last_prefix = buffer + last_prefix ;

	        if(verbose)
	        	std::cout << "Insert: " << last_prefix << " : " << suffixient_set[i] << std::endl; 
	        insert_prefix(last_prefix,suffixient_set[i]);
	        tot_inserted_char += last_prefix.size();
	        tot_inserted_keywords++;
	    }

	    file.close();
	    if(verbose)
	    {
	    	//print_trie();
	    	std::cout << "Tot inserted symbols: " << tot_inserted_char << std::endl;
	    	std::cout << "Tot inserted keywords: " << tot_inserted_keywords << std::endl;
	    }

	    G.build(text,8);

	    return std::make_pair(tot_inserted_char, tot_inserted_keywords);
	}

	// build suffixient index by indexing the supermaximal extensions
	std::pair<safe_t,safe_t> build(const std::string &text, const std::string &suffixient,
		                         		   const std::string &lcs, bool verbose = false)
	{
		std::cout << "Constructing LCP/LCS data structure for " << text << std::endl;
	    G.build(text);
		// BUILD PREFIX-SEARCH DATA STRUCTURE //
		std::vector<std::pair<uint_t,uint_t>> keywords;
		safe_t tot_inserted_char = 0, tot_inserted_keywords = 0;

	    // Open the input files in binary mode
	    std::ifstream file1(suffixient, std::ios::binary);
	    std::ifstream file2(lcs, std::ios::binary);
	    if (!file1.is_open() or !file2.is_open()) {
	        std::cerr << "Error: Could not open input files... " << std::endl;
	        exit(1);
	    }

	    safe_t pos, lcs_val; int TEMP = 0;
	    std::vector<uint8_t> suff_buffer(5,0), lcs_buffer(5,0);
	    //while (file1.read(reinterpret_cast<char*>(&pos), 5) and
	    //		  	  file2.read(reinterpret_cast<char*>(&lcs_val), 5))
	    while (file1.read(reinterpret_cast<char*>(&suff_buffer[0]), 5) and
	    	   file2.read(reinterpret_cast<char*>(&lcs_buffer[0]), 5)
	    	  )
	    {
	    	pos = get_5bytes_uint(&suff_buffer[0]);
	    	lcs_val = get_5bytes_uint(&lcs_buffer[0]);

        lcs_val += 1;
        keywords.push_back(std::make_pair(pos-lcs_val+1,lcs_val));
	    }
	    if (not (file1.eof() or file2.eof()) )
	    	{ std::cerr << "An error occurred while reading the files." << std::endl; }
	    file1.close(); file2.close();

	    //std::sort(keywords.begin(), keywords.end());

	    file1 = std::ifstream(text, std::ios::binary);
	    if (!file1.is_open()) {
	        std::cerr << "Error: Could not open file " << text << std::endl;
	        exit(1);
	    }

	    // compute file size
	    file1.seekg(0, std::ios::end);
	    safe_t file_size = file1.tellg();
	    file1.seekg(0, std::ios::beg);
	    // initialize text buffer
	    safe_t buffer_size = std::min(safe_t(DEF_BUFFER_SIZE),file_size);
	    safe_t curr = 0, last = buffer_size;
	    std::string buffer(buffer_size,'0');
	    file1.read(&buffer[0], buffer_size);

	    for(safe_t i=0;i<keywords.size();++i)
	    {
	    	safe_t start = keywords[i].first;
	    	safe_t len   = keywords[i].second;

	    	if(not contained(curr,last,start,start+len))
	    	{
	    		file1.clear();
	    		curr = start;
	    		file1.seekg(curr, std::ios::beg);
	    		buffer_size = std::max(safe_t(DEF_BUFFER_SIZE),len);
	    		buffer.clear();
	    		buffer.resize(buffer_size);
	    		file1.read(&buffer[0], buffer_size);
	    		last = curr + buffer_size;
	    	}

    		std::string keyword(keywords[i].second, '0');

    		for(safe_t i=0;i<keyword.size();++i)
    			keyword[i] = buffer[start+len-i-curr-1];

      	insert_prefix(keyword,start+len-1);
      	tot_inserted_char += keyword.size();
      	tot_inserted_keywords++;
	    }
	    file1.close();

	    return std::make_pair(tot_inserted_char, tot_inserted_keywords);
	}

	// function to insert a prefix - text position pair in the z-fast trie
	template<typename U>
	inline void insert_prefix(const std::string& prefix, const U& position)
	{
		Z.insert(&prefix,position);
		
		return;
	}

	// print z-fast trie
	void print_trie(){ Z.print(); }

	/* STORE AND LOAD NOT YET FINISHED
	Ulong store(std::string outputFile)
	{
		std::ofstream out(outputFile, std::ios::binary);
		Ulong size = Z.store(out);
		out.close();

		return size;
	}
	void load(std::string inputFile)
	{
		std::ifstream in(inputFile, std::ios::binary);
		Z.load(in);
		in.close();

		return ;
	} */

	Ulong locate_prefix(std::string pattern){ return Z.locatePrefix(pattern); }

	std::tuple<Ulong,Int,bool> locate_longest_prefix(std::string pattern){
		return Z.locateLongestPrefix(pattern);
	}

	std::string get_longest_match_seq(std::string pattern){
		return Z.getLongestMatchSeq(pattern);
	}

	void find_MEMs(string& pattern, std::ofstream& output)
	{
		size_t i = 0, l = 0, m = pattern.size();
		int64_t pstart = 0;
		while(i < m)
		{
			std::string right_max_substr = pattern.substr(pstart,(i+1)-pstart);
			std::reverse(right_max_substr.begin(),right_max_substr.end());
			auto j = locate_longest_prefix(right_max_substr);
			size_t b = G.LCS(pattern,i,std::get<0>(j));
			if(b <= l)
			{
				output << "(" << i-l << "," << l << ") ";
				pstart = i-l+1;
			}
			size_t f = G.LCP(pattern,i+1,std::get<0>(j)+1);
			i = i + f + 1;
			l = b + f;
		}
		output << "(" << i-l << "," << l << ")" << std::endl;
	}

	void find_MEMs_v2(string& pattern, std::ofstream& output)
	{
		size_t i = 0, l = 0, m = pattern.size();
		int64_t pstart = 0;
		while(i < m)
		{
			std::string right_max_substr = pattern.substr(pstart,(i+1)-pstart);
			std::reverse(right_max_substr.begin(),right_max_substr.end());
			auto j = locate_longest_prefix(right_max_substr);
			size_t b = std::get<1>(j);
			if(not std::get<2>(j)){ b += G.LCS(pattern,i-b,std::get<0>(j)-b); }
			if(b <= l)
			{
				output << "(" << i-l << "," << l << ") ";
				pstart = i-l+1;
			}
			size_t f = G.LCP(pattern,i+1,std::get<0>(j)+1);
			i = i + f + 1;
			l = b + f;
		}
		output << "(" << i-l << "," << l << ")" << std::endl;
	}

	void locate_fasta(std::string patternFile)
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
				find_MEMs_v2(line,output);
				c += line.size();
			}
			else{ output << line << std::endl; }
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

private:
	// prefix search data structure
	prefixSearch Z;
	// lcs/lcp data structure
	lcpDS G; // TO INCLUDE SLP index
};

}

#endif