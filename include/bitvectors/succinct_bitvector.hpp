/*
 * Construction of siccinct bitvectors supporting select queries
 * 
 * This code is taken from https://github.com/nicolaprezza/r-index/blob/master/internal/succinct_bit_vector.hpp
 *
 */

#ifndef SUCCINCT_BITVECTOR_HPP_
#define SUCCINCT_BITVECTOR_HPP_

#include <cassert>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <common_.hpp>

namespace suffixient{

class succinct_bitvector{

public:

	/*
	 * empty constructor. Initialize bitvector with length 0.
	 */
	succinct_bitvector(){}

	/*
	 * constructor. build bitvector given a vector of bools
	 */
	succinct_bitvector(std::vector<bool_t>& b)
	{
		bv = sdsl::bit_vector(b.size());

		for(usafe_t i=0;i<b.size();++i){ bv[i] = b[i]; }

		select1 = sdsl::bit_vector::select_1_type(&bv);
		// rank1   = sdsl::bit_vector::rank_1_type(&bv);
	}

	succinct_bitvector & operator= (const succinct_bitvector & other)
	{
		bv = sdsl::bit_vector(other.bv);
		// rank1 = sdsl::bit_vector::rank_1_type(&bv);
		select1 = sdsl::bit_vector::select_1_type(&bv);

	    return *this;
	}

	/*
	 * argument: position i in the bitvector
	 * returns: bit in position i
	 * only access! the bitvector is static.
	 */
	bool_t operator[](uint_t i) const
	{
		assert(i<size());
		return bv[i];
	}

	/*
	 * argument: position i in the bitvector
	 * returns: number of bits equal to 1 before position i excluded
	 *//*
	uint_t rank(uint_t i)
	{
		assert(i<=size());
		return rank1(i);
	}*/

	/*
	 * argument: integer i>=0
	 * returns: position of the i-th one in the bitvector. i starts from 0!
	 */
	uint_t select(uint_t i) const
	{
		assert(i<number_of_1());
		return select1(i+1);//in sd_vector, i starts from 1
	}

	/*
	* returns: size of the bitvector
	*/
	uint_t size() const { return bv.size(); }

	/*
	 * returns: number of 1s in the bitvector
	 */
	// uint_t number_of_1(){ return rank1(size()); }

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out)
	{
		uint_t size=0;

		size += bv.serialize(out);

		assert(size>0);

		return size;
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in)
	{
		bv.load(in);
		// rank1 = sdsl::bit_vector::rank_1_type(&bv);
		select1 = sdsl::bit_vector::select_1_type(&bv);
	}

private:

	// compressed bit vector
  	sdsl::bit_vector bv;
  	// rank data structure
  	// sdsl::bit_vector::rank_1_type rank1;
  	// select data structure
  	sdsl::bit_vector::select_1_type select1;
};

}

#endif // SUCCINCT_BITVECTOR_HPP_