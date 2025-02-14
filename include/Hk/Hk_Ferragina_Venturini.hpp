#ifndef __HK_FERRAGINA_VENTURINI__HPP__
#define __HK_FERRAGINA_VENTURINI__HPP__

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/io.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

#define HK_FV_HEADER (0x0d6c0000 + 0x0001)

template <class BV_CODE=sdsl::bit_vector, class BV_SEPARATOR=sdsl::sd_vector<>, class RANK_1_TYPE = typename BV_SEPARATOR::rank_1_type, class SELECT_1_TYPE = typename BV_SEPARATOR::select_1_type>
struct Hk_Ferragina_Venturini {
    void construct( const std::string& input_filename, size_t num_chars_in_block ) {
        sdsl::int_vector<8> T;
        sdsl::load_vector_from_file( T, input_filename, 1 );

        N = T.size();
        CharsInBlock = num_chars_in_block;
        int num_blocks = (N+num_chars_in_block-1)/num_chars_in_block;

        if( N % num_chars_in_block != 0 ) {
            T.resize( num_blocks * num_chars_in_block );
        }

        std::vector<size_t> ID( num_blocks );
        for( size_t i = 0; i < num_blocks; ++i ) {
            ID[i] = i;
        }

        for( size_t k = 0; k < num_chars_in_block; ++k ) {
            std::vector<size_t> count( 256, 0 );
            for( size_t i = 0; i < num_blocks; ++i ) {
                ++count[ T[ ID[i]*num_chars_in_block + k ] ];
            }
            std::vector<size_t> offset( 256 );
            offset[0] = 0;
            for( size_t i = 1; i < 256; ++i ) {
                offset[i] = offset[i-1] + count[i-1];
            }
            std::vector<size_t> _ID( num_blocks );
            for( size_t i = 0; i < num_blocks; ++i ) {
                unsigned char c = T[ ID[i]*num_chars_in_block + k ];
                _ID[ offset[c]++ ] = ID[i];
            }
            ID.swap( _ID );
        }

        std::vector<size_t> DistinctBlock;
        std::vector<size_t> Freq;
        size_t num_distinct_blocks = 1;
        size_t current_freq = 1;
        DistinctBlock.push_back( 0 );
        for( size_t i = 1; i < num_blocks; ++i ) {
            bool diff = false;
            for( size_t k = 0; k < num_chars_in_block; ++k ) {
                if(    T[ ID[i-1]*num_chars_in_block + k ] 
                    != T[ ID[i  ]*num_chars_in_block + k ] ) {

                    diff = true; break;
                }
            }
            if( diff ) {
                ++ num_distinct_blocks;
                DistinctBlock.push_back( i );
                Freq.push_back( current_freq );
                current_freq = 0;
            }
            ++ current_freq;
        } 
        Freq.push_back( current_freq );

        DistinctBlock.push_back( num_blocks );
        // NOTE: ID[Distinct[0]] is the lexicographically smallest block

        std::vector<size_t> SortedByFreq( Freq.size() );
        for( size_t i = 0; i < Freq.size(); ++i ) {
            SortedByFreq[i] = i;
        }
        std::sort( SortedByFreq.begin(), SortedByFreq.end(), 
            [&] ( size_t i, size_t j ) {
                if( Freq[i] > Freq[j] ) return true;
                if( Freq[i] < Freq[j] ) return false;
                return i < j;
            } );
        // NOTE: ID[DistinctBlock[SortedByFreq[0]]] is the most frequent block.

        std::vector<size_t> IDinv( num_blocks );
        for( size_t i = 0; i < SortedByFreq.size(); ++i ) {
            size_t b_beg = DistinctBlock[SortedByFreq[i]  ];
            size_t b_end = DistinctBlock[SortedByFreq[i]+1];

            for( size_t b = b_beg; b < b_end; ++b ) {
                IDinv[ ID[b] ] = i;
            }
        }

        // *********
        // Construct
        // *********

        size_t total_codelen = 0;
        std::vector<size_t> pos;
        for( size_t i = 0; i < IDinv.size(); ++i ) {
            size_t j = IDinv[i];
            pos.push_back( total_codelen );
            total_codelen += sdsl::bits::hi(j+2);
        }

        // code
        encoded.resize( total_codelen );
        size_t code_offset = 0;
        for( size_t i = 0; i < IDinv.size(); ++i ) {
            size_t j = IDinv[i];
            code_offset += write_code(code_offset,j);
        }

        // separator
        sdsl::bit_vector b(total_codelen+1, 0);
        for( size_t i = 0; i < pos.size(); ++i ) {
            b[pos[i]] = 1;
        }
        b[total_codelen] = 1;

        bv = BV_SEPARATOR(b);
        rank = RANK_1_TYPE(&bv);
        select = SELECT_1_TYPE(&bv);

        // table
        table.resize( num_chars_in_block * Freq.size() ); // Freq.size(): number of distinct blocks
        size_t table_write_offset = 0;
        for( size_t i = 0; i < Freq.size(); ++i ) {
            size_t block_id = ID[DistinctBlock[SortedByFreq[i]]];
            for( size_t j = 0; j < num_chars_in_block; ++j ) {
                table[ table_write_offset++ ] = T[ block_id * num_chars_in_block + j ];
            }
        }
    }

    unsigned char extract( size_t i ) const { 
        if( i >= N ) return '\0';
        unsigned int b = fetch_block_containing_pos_(i);
        return table[b*CharsInBlock+(i%CharsInBlock)];
    }

    std::string extract( size_t i, size_t j ) const { 
    // extract T[i..j]
        std::string ret;
        if( j >= N ) j = N-1;
        if( i > j ) return ret;
        ret.reserve( j-i );
        while( i <= j ) {
            ret.push_back( extract(i++) );
        }
        return ret;
    }
    
    size_t LCP( const std::string& P, size_t p_i, size_t t_j ) {
    // returns lcp( P[p_i..|P|-1], T[t_j..|T|-1] )
        int l = 0;
        const char *p = P.c_str() + p_i;
        size_t b = fetch_block_containing_pos_(t_j);
        size_t remaining_char = CharsInBlock-(t_j%CharsInBlock);

        while(1) {
            const size_t table_i = b*CharsInBlock+(t_j%CharsInBlock);
            if( t_j >= N || *(p+l) != table[table_i] ) {
                break;
            }
            t_j++;
            l++;
            if( !--remaining_char ) {
                b = fetch_block_containing_pos_(t_j);
                remaining_char = CharsInBlock;
            }
        }
        return l;
    }

    size_t LCS( const std::string& P, size_t p_i, size_t t_j ) {
    // returns lcs( P[0..p_i], T[0..t_j] )
        size_t l = 0;
        const char *p = P.c_str() + p_i;
        size_t b = fetch_block_containing_pos_(t_j);
        size_t remaining_char = (t_j%CharsInBlock)+1;

        while(1) {
            const size_t table_i = b*CharsInBlock+(t_j%CharsInBlock);
            if( t_j >= N || l > p_i || *(p-l) != table[table_i] ) {
                break;
            }
            l++;
            if( t_j-- == 0 ) break;
            if( !--remaining_char ) {
                b = fetch_block_containing_pos_(t_j);
                remaining_char = CharsInBlock;
            }
        }
        return l;
    }

    std::pair<size_t,unsigned char> LCS_char( const std::string& P, size_t p_i, size_t t_j ) {
    // returns ( lcs( P[0..p_i], T[0..t_j] ), first mismatch char on T )
        size_t l = 0;
        const char *p = P.c_str() + p_i;
        size_t b = fetch_block_containing_pos_(t_j);
        size_t remaining_char = (t_j%CharsInBlock)+1;

        while(1) {
            if( t_j >= N || l > p_i ) {
                return std::make_pair( l, (unsigned char)-1 );
            }
            const size_t table_i = b*CharsInBlock+(t_j%CharsInBlock);
            if( *(p-l) != table[table_i] ) {
                return std::make_pair( l, (unsigned char)table[table_i] );
            }
            ++l;
            if( t_j-- == 0 ) {
                return std::make_pair( l, (unsigned char)-1 );
            }
            if( !--remaining_char ) {
                b = fetch_block_containing_pos_(t_j);
                remaining_char = CharsInBlock;
            }
        }
    }

    void build( const std::string& input_filename, size_t num_chars_in_block = 4 ) {
        construct( input_filename, num_chars_in_block );
        std::string index_file = input_filename + ".hkfv";
        std::ofstream fout( index_file, std::ios::binary );
        this->serialize( fout );
        fout.close();
    }

    size_t serialize( std::ostream& out ) const {
        size_t ret = 0;
        size_t header = HK_FV_HEADER;
        ret += sdsl::serialize(header,out);
        ret += sdsl::serialize(N,out);
        ret += sdsl::serialize(CharsInBlock,out);
        ret += sdsl::serialize(table,out);
        ret += sdsl::serialize(encoded,out);
        ret += sdsl::serialize(bv,out);
        ret += sdsl::serialize(rank,out);
        ret += sdsl::serialize(select,out);
        return ret;
    }

    size_t size( void ) const {
        return sizeof(size_t) // (= header-size)
             //// size of each component:
             + sizeof(N)
             + sizeof(CharsInBlock)
             + sdsl::size_in_bytes(table)
             + sdsl::size_in_bytes(encoded)
             + sdsl::size_in_bytes(bv)
             + sdsl::size_in_bytes(rank)
             + sdsl::size_in_bytes(select);
    }

    bool load( std::istream& in ) {
        size_t header;
        sdsl::read_member(header,in);
        if( header != HK_FV_HEADER ) return false;
        sdsl::read_member(N,in);
        sdsl::read_member(CharsInBlock,in);
        table.load(in);
        encoded.load(in);
        bv.load(in);
        rank.load(in,&bv);
        select.load(in,&bv);
        return true;
    }

    bool load( const std::string& index_filename ) {
        std::ifstream fin( index_filename, std::ios::binary );
        if( !fin ) return false;
        bool ret = load( fin );
        fin.close();
        return ret;
    }

    size_t N;
    size_t CharsInBlock;
    sdsl::int_vector<8> table;
    BV_CODE encoded;
    BV_SEPARATOR bv;
    RANK_1_TYPE rank;
    SELECT_1_TYPE select;

    private:

    size_t write_code( size_t offset, size_t x ) {
        int L = sdsl::bits::hi(x+=2);
        int l = L;
        while( l > 0 ) {
            encoded[offset++] = (x&1);
            x>>=1;
            --l;
        }
        return L;
    }

    size_t fetch_block_containing_pos_( size_t i ) const {
        size_t b_i = i/CharsInBlock;
        size_t e_i = select(b_i+1);
        size_t e_j = select(b_i+2);
        size_t v = ((1<<(e_j-e_i)) | encoded.get_int( e_i, e_j-e_i ))-2;  
        return v;
    }
};

#endif
