/* static_compressed_text.cpp
 *  : modified from static_selfindex.cpp
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

#include "static_compressed_text.h"
#include <cstring>

namespace lz77index{

static_compressed_text::static_compressed_text() {
	tlen = 0;
    sigma_mapper = NULL;
}
static_compressed_text::~static_compressed_text(){
    if(sigma_mapper!=NULL)delete sigma_mapper;
}
unsigned int static_compressed_text::length() {
    return tlen;
}
unsigned char* static_compressed_text::display(unsigned int start, unsigned int end){
    if(end >= this->tlen-1){//we don't want to extract the trailing '\0'
        end = this->tlen-2;
    }
    if(start>end)return NULL;
    unsigned char* sstr = _display(start,end);
    if(sstr!=NULL){
        //std::cout<<sstr<<std::endl;
        if(sigma_mapper->unmap(sstr,end-start+1)!=0){
            delete [] sstr;
            return NULL;
        }
    }
    return sstr;
}
static_compressed_text* static_compressed_text::load(FILE *fp) {
    unsigned char rd;
    if(fread(&rd,sizeof(unsigned char),1,fp)!=1) return NULL;
    fseek(fp,-1*sizeof(unsigned char),SEEK_CUR);
    switch(rd) {
        case COMPRESSED_TEXT_LZ77_HDR: return static_compressed_text_lz77::load(fp);
        //case COMPRESSED_TEXT_NONE_HDR: return static_compressed_text_none::load(fp);
        //case COMPRESSED_TEXT_LZEND_HDR: return static_compressed_text_lzend::load(fp);
        default:;
    }
    return NULL;
}
static_compressed_text* static_compressed_text::load(const char* filename) {
    FILE* fd = fopen(filename,"r");
    static_compressed_text* ret = load(fd);
    fclose(fd);
    return ret;
}

}

