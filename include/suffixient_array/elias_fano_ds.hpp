/*
 * Construction of the Elias-Fano compressed bitvectors
 * 
 * This code is taken from https://github.com/davidecenzato/extended_r-index/blob/main/sd_vector.hpp
 *
 */

#ifndef ELIAS_FANO_DS_HPP_
#define ELIAS_FANO_DS_HPP_

#include <cassert>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <common.hpp>

namespace suffixient{

class elias_fano_ds{

public:

	elias_fano_ds(){}

	elias_fano_ds(std::vector<uint_t>& onset, uint_t bsize)
	{
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset.size());
		for(auto idx: onset){ builder.set(idx); }
		bv = sdsl::sd_vector<>(builder);
		onset.clear(); 
		// set bitvector len
		u = bv.size();
	}
	// 2nd constructor
	elias_fano_ds(std::ifstream& onset, uint_t bsize)
	{
		// compute onset size
		onset.seekg(0, std::ios::end);
		size_t onset_size = onset.tellg()/sizeof(uint_t);
    	onset.seekg(0, std::ios::beg);
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset_size);
		// compute bitvector builder
		uint_t currBitPos = 0;
		for(uint_t i=0;i<onset_size;++i){
			// get new onset pos
			onset.read(reinterpret_cast<char*>(&currBitPos), sizeof(uint_t));
			builder.set(currBitPos);
		}
		bv = sdsl::sd_vector<>(builder);
		// close stream
		onset.close(); 
		// set bitvector len
		u = bv.size();
	}

	// 3nd constructor
	elias_fano_ds(std::ifstream& onset, uint_t bsize, int isize = 4)
	{
		// compute onset size
		onset.seekg(0, std::ios::end);
		size_t onset_size = onset.tellg() / isize;
    	onset.seekg(0, std::ios::beg);
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset_size);
		// compute bitvector builder
		uint64_t currBitPos = 0;
		for(uint_t i=0;i<onset_size;++i){
			// get new onset pos
			onset.read(reinterpret_cast<char*>(&currBitPos), isize);
			builder.set((uint_t)currBitPos);
		}
		bv = sdsl::sd_vector<>(builder);
		// close stream
		onset.close(); 
		// set bitvector len
		u = bv.size();
	}

	void construct_rank_ds()
	{
		assert(bv.size() > 0);
		sdsl::util::init_support(rank1_,&bv);
	}

	void construct_select_ds()
	{
		assert(bv.size() > 0);
		sdsl::util::init_support(select1_,&bv);
	}

	uint_t size(){ return u; }
	
	uint_t rank1(uint_t i){ return rank1_(i); }

	uint_t select1(uint_t i){ return select1_(i+1); }

	uint_t at(uint_t i){ return bv[i]; }

	uint_t gapAt(uint_t i)
	{
		if(i==0){ return select1(0)+1; }
		return select1(i)-select1(i-1);
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in)
	{
		in.read((char*)&u, sizeof(u));
		bv.load(in);

		construct_rank_ds();
		construct_select_ds();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out)
	{
		uint_t w_bytes = 0;

		out.write((char*)&u, sizeof(u));

		w_bytes += sizeof(u);

		if(u==0) return w_bytes;

		w_bytes += bv.serialize(out);

		return w_bytes;
	}

private:
	// bitvector lev
	uint_t u = 0;
	// compressed bit vector
  	sdsl::sd_vector<> bv;
  	// rank data structure
  	sdsl::rank_support_sd<> rank1_;
  	// select data structure
  	sdsl::select_support_sd<> select1_;
};

}

#endif // ELIAS_FANO_DS_HPP_