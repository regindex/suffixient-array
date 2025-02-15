/* static_compressed_text_lz77,h
 * : modified from selfindex_lz77.h
 * Copyright (C) 2009, Sebastian Kreft C., all rights reserved.
 *
 * static_compressed_text_lz77 definition
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef __STATIC_COMPRESSED_TEXT_LZ77_H__
#define __STATIC_COMPRESSED_TEXT_LZ77_H__

#include <cstring>
#include <vector>
#include "static_compressed_text.h"
#include "deltacodes.h"
#include "dfuds.h"
#include "static_range.h"
#include "static_permutation.h"
#include "static_sequence_wt.h"
#include "directcodes.h"

namespace lz77index{
/** "Selfindex" that uses BMH to find occurrences.
 *  It's not really a selfindex, it is used just for testing.
 *
 *  @author Sebastian Kreft C.
 */
class static_compressed_text_lz77 : public static_compressed_text {
    public:
        ~static_compressed_text_lz77();
        unsigned int size();
        /** Saves the index to a file */
        unsigned int save(FILE *fp);
        /** Loads the index from a file */
        static static_compressed_text_lz77* load(FILE * fp);
        static static_compressed_text_lz77* build(char* filename, char* filenameidx, unsigned char binaryRev, unsigned char binarySst, unsigned char store_sst_ids);
        unsigned char extract( unsigned i );

        size_t LCP( const std::string& P, size_t p_i, size_t t_j );
        size_t LCS( const std::string& P, size_t p_i, size_t t_j );
        std::pair<size_t,unsigned char> LCS_char( const std::string& P, size_t p_i, size_t t_j );

    protected:
        /** Returns the substring of the text in a new allocated array*/
        unsigned char* _display(unsigned int start, unsigned int end);
        /** Returns the size of the index */
        static_compressed_text_lz77();
        /** Extracts a substring from a LZ77 parse*/        
        void _charextractLZ77(unsigned char* answer, unsigned start, unsigned end, unsigned int offset_start, unsigned int offset_end);
        /** returns the start position of the source of phrase i, i starting from 0 */
        int _source(unsigned int i);        
        /** returns the length of a source*/
        unsigned int _source_length(unsigned int i);

        /** additionally implemented member functions **/
        size_t m_find_mismatch_forward( const std::string& P, size_t p_i, size_t start, size_t end );
        std::pair<size_t, unsigned char> m_find_matchspan_backward( const std::string& P, size_t p_i, size_t start, size_t end );

        /** instance variables*/
        unsigned int parsing_len;
        unsigned int* trailing_char;
        unsigned char trailing_char_bits;
        DeltaCodes* phrases;
        DeltaCodes* sources;
        basics::static_permutation* perm;
        static_sequence_wt* wt_depth;
        //configuration variables
        unsigned char binaryRev;
        unsigned char binarySst;
        unsigned char store_sst_ids;
};

}
#endif /* _STATIC_SELFINDEX_LZ77_H */
