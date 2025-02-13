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

class elias_fano_bitvector{

public:

	elias_fano_bitvector(){}

	elias_fano_bitvector(std::vector<uint_t>& onset, uint_t bsize)
	{
		m = onset.size();
		// construct the compressed bitvector
		sdsl::sd_vector_builder builder(bsize,onset.size());
		for(auto idx: onset){ builder.set(idx); }
		bv = sdsl::sd_vector<>(builder);
		// onset.clear(); 
		// set bitvector len
		n = bv.size();
	}
	// 2nd constructor
	elias_fano_bitvector(std::ifstream& onset, uint_t bsize)
	{
		// compute onset size
		onset.seekg(0, std::ios::end);
		size_t onset_size = m = onset.tellg()/sizeof(uint_t);
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
		n = bv.size();
	}

	// 3nd constructor
	elias_fano_bitvector(std::ifstream& onset, uint_t bsize, int isize = 4)
	{
		// compute onset size
		onset.seekg(0, std::ios::end);
		size_t onset_size = m = onset.tellg() / isize;
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
		n = bv.size();
	}

	elias_fano_bitvector & operator= (const elias_fano_bitvector & other)
	{
		n = other.n;
		m = other.m;
		bv = sdsl::sd_vector<>(other.bv);
		rank1_ = sdsl::rank_support_sd<>(other.rank1_);
		select1_ = sdsl::select_support_sd<>(other.select1_);

	    return *this;
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

	uint_t size() const { return n; }

	uint_t no_ones() const { return m; }
	
	uint_t rank1(uint_t i) const { return rank1_(i); }

	uint_t select1(uint_t i) const { return select1_(i+1); }

	uint_t at(uint_t i) const { return bv[i]; }

	uint_t gapAt(uint_t i) const
	{
		if(i==0){ return select1(0)+1; }
		return select1(i)-select1(i-1);
	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in)
	{
		in.read((char*)&n, sizeof(n));
		in.read((char*)&m, sizeof(m));
		bv.load(in);

		construct_rank_ds();
		construct_select_ds();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	uint_t serialize(std::ostream& out) const
	{
		uint_t w_bytes = 0;

		out.write((char*)&n, sizeof(n));
		out.write((char*)&m, sizeof(m));

		w_bytes += sizeof(n) + sizeof(m);

		if(n == 0) return w_bytes;

		w_bytes += bv.serialize(out);

		return w_bytes;
	}

private:
	// bitvector length
	uint_t n = 0;
	// no. bits set
	uint_t m = 0;
	// compressed bit vector
  	sdsl::sd_vector<> bv;
  	// rank data structure
  	sdsl::rank_support_sd<> rank1_;
  	// select data structure
  	sdsl::select_support_sd<> select1_;
};

}

#endif // ELIAS_FANO_DS_HPP_