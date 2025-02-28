#ifndef __RLZ_DNA_HPP__
#define __RLZ_DNA_HPP__

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sdsl/construct.hpp>
#include <cassert>

#define RLZ_HEADER (0x0e8f0000 + 0x0001)

template < class SD_VECTOR = sdsl::sd_vector<> >
struct RLZ_DNA { 
    typedef typename SD_VECTOR::rank_1_type   rank_t;
    typedef typename SD_VECTOR::select_1_type select_t;

    struct builder {
        typedef std::tuple<size_t,size_t,unsigned char> phrase_t;

        builder( const sdsl::int_vector<8>* _text )
        : text(*_text), prefix_len(0) {
            ;
        }

        builder( const sdsl::int_vector<8>* _text, size_t _prefix_len )
        : text(*_text), prefix_len(_prefix_len) {
            if( text.size() < prefix_len ) prefix_len = text.size();
            sdsl::int_vector<8> text_prefix; text_prefix.resize( prefix_len );
            for( size_t i = 0; i < prefix_len; ++i ) {
                text_prefix[i] = text[i];
            }
            sdsl::append_zero_symbol(text_prefix);
            sdsl::algorithm::calculate_sa<>( (const unsigned char*) text_prefix.data(), text_prefix.size(), sa );

            phrase = do_parse( prefix_len, prefix_len );
        }

        builder& operator = ( const builder& _b ) {
            assert( &text == &_b.text );
            sa = _b.sa;
            prefix_len = _b.prefix_len;
            phrase = _b.phrase;
            return *this;
        }

        private:
        size_t refine_lower_bound( size_t s, size_t e, size_t j, unsigned char x ) {
            size_t n = sa[0];
            while( s < e ) {
                size_t m = (s+e)/2;
                if( sa[m]+j >= n || text[ sa[m]+j ] < x ) {
                    s=m+1;
                } else {
                    e=m;
                }
            }
            assert( s == e );
            return s;
        }

        size_t refine_upper_bound( size_t s, size_t e, size_t j, unsigned char x ) {
            size_t n = sa[0];
            while( s < e ) {
                size_t m = (s+e)/2;
                if( sa[m]+j >= n || text[ sa[m]+j ] <= x ) {
                    s=m+1;
                } else {
                    e=m;
                }
            }
            assert( s == e );
            return s;
        }

        std::pair<size_t,size_t> match( size_t j0 = 0, size_t max_phrase_length = 0 ) {
            size_t n = sa[0];
            size_t m = text.size();
            if( max_phrase_length != 0 && j0 + max_phrase_length < m ) m = j0 + max_phrase_length;
            size_t s = 0;
            size_t e = sa.size();
            size_t max_lcp = 0;
            size_t max_i = text[j0];
            size_t j = 0;
            while(j0+j<m) {
                s = refine_lower_bound( s, e, j, text[j0+j] );
                e = refine_upper_bound( s, e, j, text[j0+j] );
                if( s < e ) { ++j; max_i = sa[s]; max_lcp++; }
                else break;
            }
            return std::make_pair( max_i, max_lcp );
        }

        std::vector< std::tuple<size_t,size_t,unsigned char> > do_parse( size_t j, size_t max_phrase_length ) {
            size_t m = text.size();
            std::vector< std::tuple<size_t,size_t,unsigned char> > ret;    
            while( j < m ) {
                std::pair<size_t,size_t> r = match( j, max_phrase_length );
                j += r.second;
                if( j < m ) ret.push_back( std::make_tuple( r.first, r.second, (unsigned char)text[j] ) );
                else        ret.push_back( std::make_tuple( r.first, r.second, (unsigned char)0    ) );
                ++j;
            }
            return ret;
        }

        const sdsl::int_vector<8>& text;
        sdsl::int_vector<32> sa;
        size_t prefix_len;
        std::vector< phrase_t > phrase;

        public:
        size_t total_length    ( void ) const { return text.size(); }
        size_t reference_length( void ) const { return prefix_len; }
        size_t num_phrases     ( void ) const { return phrase.size(); }
        size_t phrase_offset   ( size_t i ) const { return std::get<0>(phrase[i]); }
        size_t phrase_length   ( size_t i ) const { return std::get<1>(phrase[i]); }
        size_t phrase_char     ( size_t i ) const { return std::get<2>(phrase[i]); }
    };

    struct bit_packed_DNA_string {
        static uint64_t pack_uint64( uint64_t x ) {
            x = (((x>>1)^x)>>1)&0x0303030303030303ull;
            x |= x>>6;
            x |= x>>12;
            x = ((x>>24)&0xff00) | (x&0xff);
            return x;
        }

        static uint16_t get_packed_chunk_from_vec_8_unsafe( const sdsl::int_vector<8>& x, size_t i ) {
            return pack_uint64( *(uint64_t*)x.data()+i );
        }

        static uint16_t get_packed_chunk_from_vec_8( const sdsl::int_vector<8>& x, size_t i ) {
            if( i+8 <= x.size() ) { return pack_uint64( *(uint64_t*)((char*)x.data()+i) ); }
            else {
                uint64_t v = 0;
                uint64_t o = 0;
                while( i < x.size() ) {
                    v |= (((uint64_t)x[i++]) << ((o++)<<3) );
                }
                return pack_uint64( v );
            }
        }

        void build( const sdsl::int_vector<8>& t, size_t i, size_t j ) {
            len = j-i;
            seq.resize( (j-i+31)/32 );
            size_t x   = 0;
            size_t cnt = 0;
            size_t write_i = 0;
            while( i < j ) {
                uint64_t chunk = get_packed_chunk_from_vec_8(t,i);
                x |= chunk << (cnt<<4);
                if( ++cnt == 4 ) {
                    seq[write_i++] = x;
                    x = cnt = 0;
                }
                i += 8;
            }
            if( cnt ) {
                seq[write_i] = x;
            }
        }

        char extract( size_t i ) const {
            if( i >= len ) return '\0';
            return ( (0x54474341u) >> 
                     ( ((seq[i>>5]>>((i&31)<<1))&3)<<3) 
                   ) & 0xff;
        }

        size_t serialize( std::ostream& out ) const {
            size_t ret = 0;
            ret += sdsl::serialize( len, out );
            ret += sdsl::serialize( seq, out );
            return ret;
        }

        size_t size( void ) const {
            size_t ret = 0;
            ret += sizeof(len);
            ret += sdsl::size_in_bytes(seq);
            return ret;
        }

        bool load( std::istream& in ) {
            sdsl::read_member( len, in );
            seq.load( in );
            return true;
        }

        bool load( const std::string& in_filename ) {
            std::ifstream fin( in_filename, std::ios::binary );
            if( !fin ) return false;
            bool ret = load( fin );
            fin.close();
            return ret;
        }

        size_t len; // number of characters
        sdsl::int_vector<64> seq;
    };

    void build( const std::string& input_filename, double epsilon = 0.1, size_t __prefix_len = 0 ) {
        sdsl::int_vector<8> text;
        sdsl::load_vector_from_file(text, input_filename, 1);

        size_t text_len   = text.size();
        size_t prefix_len = __prefix_len;

        size_t p_b, prefix_size_upper_bound;
        if( prefix_len == 0 ) {
            p_b = std::min( (((uint64_t)1)<<20), text_len/1024+1 );
            prefix_size_upper_bound = std::min( text_len, ((uint64_t)1)<<(64-8)) ;
        } else {
            p_b = prefix_len; 
            prefix_size_upper_bound = p_b+1;
        }
        size_t min_size = 0;
        builder b(&text);
        while( p_b < prefix_size_upper_bound ) { 
            // ASSUME: prefix is short enough so that (offset,char) can be encoded in 64bits

            builder b_(&text, p_b);
            size_t m = b_.num_phrases();

            size_t size_estimated = p_b * 2  // size of bit-pack string
                                  + m * ( 2 + (log(text_len) - log(m))/log(2) + 64 ); // RLZ

            if( min_size == 0 || size_estimated < min_size ) {
                min_size = size_estimated;
                b = b_;
                prefix_len = p_b;
            }
            if( size_estimated > min_size ) {
                break;
            }
            
            if( p_b == (size_t)((1+epsilon)*p_b) ) ++p_b;
            else p_b = (1+epsilon)*p_b;
        }

        total_length = text_len;
        reference.build( text, 0, prefix_len );

        parse_info.resize( b.num_phrases() );
        sdsl::sd_vector_builder sd_builder( text_len-prefix_len+2, b.num_phrases()+1 ); 
        sd_builder.set(0);
        size_t relative_offset = 0;
        for( size_t i = 0; i < b.num_phrases(); ++i ) {
            relative_offset += b.phrase_length(i)+1;
            sd_builder.set( relative_offset );
            parse_info[i] = (((uint64_t)b.phrase_offset(i)) << 8 ) | (uint8_t)b.phrase_char(i);
        }
        boundary = SD_VECTOR( sd_builder );
        b_rank   = rank_t( &boundary );
        b_select = select_t( &boundary );

        std::ofstream fout( input_filename + ".rlz", std::ios::binary );
        serialize( fout );
        fout.close();
    }

    size_t serialize( std::ostream& out ) const {
        size_t ret = 0;
        uint64_t header = RLZ_HEADER;
        ret += sdsl::serialize( header, out );
        ret += sdsl::serialize( total_length, out );
        ret += reference.serialize( out );
        ret += sdsl::serialize( boundary  , out );
        ret += sdsl::serialize( b_rank    , out );
        ret += sdsl::serialize( b_select  , out );
        ret += sdsl::serialize( parse_info, out );
        return ret;
    }

    size_t size( void ) const {
        size_t ret = 0;
        ret += sizeof( uint64_t );
        ret += sizeof( total_length );
        ret += reference.size();
        ret += sdsl::size_in_bytes( boundary     );
        ret += sdsl::size_in_bytes( b_rank       );
        ret += sdsl::size_in_bytes( b_select     );
        ret += sdsl::size_in_bytes( parse_info   );
        return ret;
    }

    unsigned char extract( size_t i ) const {
        if( i >= total_length ) return '\0';
        if( i < reference.len ) return reference.extract( i );
        i -= reference.len;
        size_t blk_id    = b_rank(i+1)-1;
        size_t p_info    = parse_info[blk_id];
        size_t offset    = p_info >>   8;
        unsigned char ch = p_info & 0xff;

        size_t curr_begin = b_select( blk_id+1 );
        size_t next_begin = b_select( blk_id+2 );

        assert( curr_begin <= i );
        assert( i <= next_begin );

        if( next_begin == i+1 ) return ch;
        return reference.extract( offset + i - curr_begin ); 
    }

    size_t LCP( const std::string& P, size_t p, size_t t ) const {
        size_t m = P.size();
        size_t n = total_length;
        size_t l = 0;
        while( p+l < m && t+l < n ) {
            unsigned char ch = extract(t+l);
            if( P[p+l] != ch ) break;
            ++l;
        }
        return l;
    }

    size_t LCS( const std::string& P, size_t p, size_t t ) const {
        size_t l = 0;
        while( l <= p && l <= t ) {
            unsigned char ch = extract(t-l);
            if( P[p-l] != ch ) break;
            ++l;
        }
        return l;
    }

    std::pair<size_t,unsigned char> LCS_char( const std::string& P, size_t p, size_t t ) const {
        size_t l = 0;
        while( l <= p && l <= t ) {
            unsigned char ch = extract(t-l);
            if( P[p-l] != ch ) return std::make_pair(l,ch);
            ++l;
        }
        return std::make_pair(l,-1);
    }

    bool load( std::ifstream& in ) {
        uint64_t header;
        sdsl::read_member( header, in );
        if( header != RLZ_HEADER ) return false;
        sdsl::read_member( total_length, in );
        reference .load( in );
        boundary  .load( in );
        b_rank    .load( in, &boundary );
        b_select  .load( in, &boundary );
        parse_info.load( in );
        return !!in;
    }

    bool load( const std::string& index_filename ) {
        std::ifstream fin( index_filename, std::ios::binary );
        if( !fin ) return false;
        bool ret = load( fin );
        fin.close();
        return ret;
    }

    uint64_t     total_length;
    bit_packed_DNA_string reference;
    SD_VECTOR    boundary;
    rank_t       b_rank;
    select_t     b_select;
    sdsl::int_vector<64>  parse_info;
};

#endif
