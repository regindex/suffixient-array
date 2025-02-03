/* static_compressed_text.h
 *  : modified from static_selfindex.h
 * Copyright (C) 2009, Sebastian Kreft C., all rights reserved.
 *
 * static_compressed_text definition
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

#ifndef __STATIC_COMPRESSED_TEXT_H__
#define __STATIC_COMPRESSED_TEXT_H__

#include <cstdio>
#include <iostream>
#include <vector>
#include "mapper.h"

#define COMPRESSED_TEXT_LZ77_HDR 2
//#define SELFINDEX_NONE_HDR 3
//#define SELFINDEX_LZEND_HDR 4

namespace lz77index{
/** Base class for selindexes
 * 
 *  @author Sebastian Kreft
 */
class static_compressed_text {
    public:
        static_compressed_text();
        virtual ~static_compressed_text();
        /** Returns the substring of the text in a new allocated array*/
        unsigned char* display(unsigned int start, unsigned int end);
        /** Returns the length of the text*/
        unsigned int length();
        /** Returns the size of the index */
        virtual unsigned int size()=0;
        /** Saves the index to a file */
        virtual unsigned int save(FILE *fp)=0;
        /** Loads the index from a file */
        static static_compressed_text* load(FILE * fp);
        static static_compressed_text* load(const char* filename);
    protected:
        virtual unsigned char* _display(unsigned int start, unsigned int end)=0;
        unsigned int tlen;
        mapper* sigma_mapper;        
};

}

#include "static_compressed_text_lz77.h"
//#include "static_compressed_text_none.h"
//#include "static_compressed_text_lz77.h"
//#include "static_compressed_text_lzend.h"

#endif /* _STATIC_SELFINDEX_H */
