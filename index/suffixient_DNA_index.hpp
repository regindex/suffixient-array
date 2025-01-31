// Copyright (c) 2024, REGZ.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  suffixient_index: DESCRIPTION
 */

#ifndef SUFFIXIENT_DNA_INDEX_HPP_
#define SUFFIXIENT_DNA_INDEX_HPP_

#include <chrono>
#include <malloc_count.h>
// suffixient array data structure
#include <suffixient_array_baseline.hpp>
// Text oracle data structure
#include <LZ77_Kreft_Navarro_index.hpp>

namespace suffixient{

template<class suffArray, class textOracle>
class suffixient_DNA_index{

public:
	// empty constructor
	suffixient_DNA_index(){}

	// build suffixient index by indexing the supermaximal extensions
	std::pair<usafe_t,usafe_t> build(const std::string &text_filepath, const uint8_t len)
	{
		// BUILD LZ77 INDEX DATA STRUCTURE //
		// std::cout << "Constructing text oracle data structure for " << text_filepath << std::endl;
	  	// O.build(text_filepath);

	  	// malloc_count_reset_peak();
	  	// std::cout << "mem peak: " << malloc_count_peak() << std::endl;
		// BUILD PREFIX-SEARCH DATA STRUCTURE //
		S.build(text_filepath, len);

	  return std::make_pair(0, 0);
	}

private:
	// suffixient array search data structure
	suffArray S;
	// text oracle data structure
	textOracle O;
};

}

#endif