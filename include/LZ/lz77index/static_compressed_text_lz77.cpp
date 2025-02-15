#include "static_compressed_text_lz77.h"

#include <algorithm>
#include <string>

#include "utils.h"
#include "utils_index.h"
#include "mapper.h"
#include "LZparser.h"
#include "LZ77.h"
#include "LZend.h"
#include "patricia.h"
#include "static_bitsequence_builder.h"
#include "directcodes.h"

static const int DELTA_SAMP = 16;
static const int PERM_EPS = 8;
static const int BITM_STEP = 20;

namespace lz77index{
static_compressed_text_lz77* static_compressed_text_lz77::build(char* filename, char* filenameidx, unsigned char binaryRev, unsigned char binarySst, unsigned char store_sst_ids){

    std::string sfile = std::string((char*)filename);
    std::cout<<"Building the index to file "<<filenameidx<<std::endl;

    FILE* lz77index_fp = fopen(filenameidx,"w");
    //0) Save the header
    unsigned char hdr = COMPRESSED_TEXT_LZ77_HDR;
    if(fwrite(&hdr,sizeof(unsigned char),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving the header"<<std::endl;
        return NULL;
    }
    //0') Save the options
    if(fwrite(&binaryRev,sizeof(unsigned char),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving option"<<std::endl;
        return NULL;
    }
    if(fwrite(&binarySst,sizeof(unsigned char),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving option"<<std::endl;
        return NULL;
    }
    if(fwrite(&store_sst_ids,sizeof(unsigned char),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving option"<<std::endl;
        return NULL;
    }
    /************************************************************************************************************************************************************************************************/
    //1) Read the file into an array
    unsigned int tlen;
    unsigned char* text = utils::readText((char*)filename,&tlen);
    if(text == NULL){
        std::cerr<<"ERROR: reading file"<<sfile<<std::endl;
        return NULL;
    }
    if(fwrite(&tlen,sizeof(unsigned int),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving the text size"<<std::endl;
        return NULL;
    }
    /************************************************************************************************************************************************************************************************/
    //2) Create the alphabet mapper
    mapper* sigma_mapper = new mapper(text,tlen);
    unsigned int sigma = sigma_mapper->sigma();
    std::cout<<"Sigma: "<<sigma<<std::endl;
    /************************************************************************************************************************************************************************************************/
    //3) Convert the text to the new alphabet and save it to a file
    sigma_mapper->map(text,tlen);
    utils::saveArray(text,tlen,(sfile+std::string(".alph")).c_str());
    /************************************************************************************************************************************************************************************************/
    //4) Save mapper
    if(sigma_mapper->save(lz77index_fp)!=0){
        std::cerr<<"ERROR: while saving the alphabet mapper"<<std::endl;
        return NULL; 
    }
    std::cout<<"(S)Mapper: "<<sigma_mapper->size()<<std::endl;
    /************************************************************************************************************************************************************************************************/
    //5) Free memory(mapper and text)
    delete sigma_mapper;sigma_mapper=NULL;
    delete[] text;
    /************************************************************************************************************************************************************************************************/
    //6) Compute LZ77(or like) parsing over the new text, it will generate three different files, .len, .start and .char
    parsing::LZparser* parser = new parsing::LZ77((sfile+std::string(".alph")).c_str());
    parser->parse();
    delete parser;
    /************************************************************************************************************************************************************************************************/    
    //7) Load .char array and save it to lz77index_fp using nlog(sigma) bits
    unsigned int parsing_len;
    unsigned char* trailing_char = (unsigned char*)utils::readArray((sfile+std::string(".alph.char")).c_str(),&parsing_len,1);
    std::cout<<"Factors: "<<parsing_len<<std::endl;

    //we will store all character up to parsing_len-1 as char[i]-1, and the last one as 0
    //because we know that always the last character is 0, thus we reduce sigma by one
    unsigned int* trailing_char_compressed = new unsigned int[basics::uint_len(parsing_len,basics::bits(sigma-1))];
    for(unsigned int i=0;i<basics::uint_len(parsing_len,basics::bits(sigma-1));i++){
        trailing_char_compressed[i] = 0;
    }
    for(unsigned int i=0;i<parsing_len-1;i++){
        basics::set_field(trailing_char_compressed, basics::bits(sigma-1), i, trailing_char[i]-1);
    }
    basics::set_field(trailing_char_compressed, basics::bits(sigma-1), parsing_len-1, 0);
    if(fwrite(&parsing_len,sizeof(unsigned int),1,lz77index_fp)!=1){
        std::cerr<<"ERROR: while saving the parsing size"<<std::endl;
        return NULL; 
    }
    if(fwrite(trailing_char_compressed,sizeof(unsigned int),basics::uint_len(parsing_len,basics::bits(sigma-1)),lz77index_fp)!=basics::uint_len(parsing_len,basics::bits(sigma-1))){
        std::cerr<<"ERROR: while saving the trailing character array"<<std::endl;
        return NULL; 
    }
    std::cout<<"(S)Trailing character: "<<sizeof(unsigned int)*basics::uint_len(parsing_len,basics::bits(sigma-1))<<std::endl;
    delete[] trailing_char;
    delete[] trailing_char_compressed;
    /************************************************************************************************************************************************************************************************/
    //8) Compute bitmap of phrases    
    unsigned int* length = (unsigned int*)lz77index::utils::readArray((sfile+".alph.len").c_str(),&parsing_len,4);
    DeltaCodes* phrases = new DeltaCodes(length,parsing_len,DELTA_SAMP);//sampling=8
    //delete [] length;
    if(phrases->save(lz77index_fp)!=0){
        std::cerr<<"ERROR: while saving the DeltaCodes bitmap"<<std::endl;
        return NULL;
    }
    std::cout<<"(S)Phrases: "<<phrases->size()<<std::endl;
    delete phrases;
    //11 Compute depth array and bitmap(it is actually a prevless structure)
    //12 compute permutation
    unsigned int len;
    length = (unsigned int*)utils::readArray((sfile+".alph.len").c_str(),&len,4);
    unsigned int* start = (unsigned int*)utils::readArray((sfile+".alph.start").c_str(),&len,4);
    for(unsigned int i=0;i<len;i++){
        length[i] -= 1;
        if(length[i]!=0){
            start[i] += 1;
        }else{
            start[i] = 0;
        }
    }
    utils_index::Phrase* aux_phrases = new utils_index::Phrase[len];
    for(unsigned int i=0;i<len;i++){
        aux_phrases[i].start_pos = start[i];
        //std::cout<<start[i]<<" ";
        aux_phrases[i].length = length[i];
        aux_phrases[i].id = i;
    }
    //std::cout<<std::endl;
    delete [] length;
    delete [] start;
    unsigned char* depth = utils_index::computeDepth(aux_phrases,len);
    /*for(unsigned int i=0;i<len;i++){
        depth[i] -= 1;
        //std::cout<<(unsigned int)depth[i]<<" ";
    }*/
    //std::cout<<std::endl;
    unsigned int sigma_d = *std::max_element(depth,depth+len)+1;
    std::cout<<"Sigma depth: "<<sigma_d<<std::endl;
    static_sequence_wt* wt_depth = new static_sequence_wt(depth,len,(unsigned char)sigma_d);
    if(wt_depth->save(lz77index_fp)!=0){
        std::cerr<<"ERROR: while saving the depth WT"<<std::endl;
        return NULL; 
    }
    std::cout<<"(S)depth: "<<wt_depth->size()<<std::endl;
    delete wt_depth;

    delete [] depth;
    unsigned int* permutation_array = utils_index::computePermutationArray(aux_phrases,len);
    basics::static_bitsequence_builder* bmb = new basics::static_bitsequence_builder_brw32(BITM_STEP);
    basics::static_permutation* perm = new basics::static_permutation_mrrr(permutation_array,len,PERM_EPS,bmb);//probar que valor funciona mejor
    delete bmb;
    if(perm->save(lz77index_fp) != 0){
        std::cerr<<"ERROR: while saving the permutation"<<std::endl;
        return NULL;
    }
    std::cout<<"(S)perm: "<<perm->size()<<std::endl;

    delete perm;
    unsigned int* bitmap_sources = utils_index::computeSourcesBitmap(aux_phrases,len);
    phrases = new DeltaCodes(bitmap_sources,len,DELTA_SAMP);//sampling=8

    if(phrases->save(lz77index_fp)!=0){
        std::cerr<<"ERROR: while saving the DeltaCodes sources bitmap"<<std::endl;
        return NULL;
    }
    std::cout<<"(S)sources: "<<phrases->size()<<std::endl;
    delete phrases;
    delete [] bitmap_sources;
    delete [] aux_phrases;
    fclose(lz77index_fp);
    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    std::cout<<"Loading the index from file "<<(filenameidx)<<std::endl;
    static_compressed_text_lz77* ret;

    lz77index_fp = fopen(filenameidx,"r");
    ret = static_compressed_text_lz77::load(lz77index_fp);
    fclose(lz77index_fp);
    return ret;
}
static_compressed_text_lz77::~static_compressed_text_lz77(){
    //if(sigma_mapper != NULL) delete sigma_mapper;
    if(trailing_char != NULL) delete [] trailing_char;
    if(phrases != NULL) delete phrases;
    if(sources != NULL) delete sources;
    if(perm != NULL)delete perm;
    if(wt_depth != NULL)delete wt_depth;
}
unsigned int static_compressed_text_lz77::size(){
    return sizeof(static_compressed_text_lz77) + 
    sigma_mapper->size() +
    sizeof(unsigned int)*basics::uint_len(parsing_len,basics::bits(sigma_mapper->sigma()-1)) +
    phrases->size() + 
    //sizeof(unsigned int)*basics::uint_len(parsing_len-1,rev_ids_bits) + 
    perm->size() + 
    sources->size();
}
/** Saves the index to a file */
unsigned int static_compressed_text_lz77::save(FILE *fp){
    return 1;
}
/** Loads the index from a file */
static_compressed_text_lz77* static_compressed_text_lz77::load(FILE * fp){
    unsigned char hdr;
    if(fread(&hdr,sizeof(unsigned char),1,fp)!=1 || hdr != COMPRESSED_TEXT_LZ77_HDR){
        std::cerr<<"ERROR: load: wrong header"<<std::endl;
        return NULL;
    }
    static_compressed_text_lz77* ret = new static_compressed_text_lz77();
    /**************************************************************************/
    /* Load Options                                                           */
    if(fread(&(ret->binaryRev),sizeof(unsigned char),1,fp)!=1){
        std::cerr<<"ERROR: load: option"<<std::endl;
        return NULL;
    }
    if(fread(&(ret->binarySst),sizeof(unsigned char),1,fp)!=1){
        std::cerr<<"ERROR: load: option"<<std::endl;
        return NULL;
    }
    if(fread(&(ret->store_sst_ids),sizeof(unsigned char),1,fp)!=1){
        std::cerr<<"ERROR: load: option"<<std::endl;
        return NULL;
    }
    /**************************************************************************/
    if(fread(&(ret->tlen),sizeof(unsigned int),1,fp) != 1){
        delete ret; return NULL;
    }    
    ret->sigma_mapper = mapper::load(fp);
    if(ret->sigma_mapper == NULL){
        std::cerr<<"ERROR: load: mapper "<<std::endl;
        delete ret; return NULL;
    }
    if(fread(&(ret->parsing_len),sizeof(unsigned int),1,fp) != 1){
        std::cerr<<"ERROR: load: parsing len "<<std::endl;
        delete ret; return NULL;
    }

    ret->trailing_char_bits = basics::bits(ret->sigma_mapper->sigma()-1);
    ret->trailing_char = new unsigned int[basics::uint_len(ret->parsing_len, ret->trailing_char_bits)];
    if(ret->trailing_char == NULL){
        std::cerr<<"ERROR: load: allocating trailing char "<<std::endl;
        delete ret; return NULL;
    }
    if(fread(ret->trailing_char,sizeof(unsigned int),basics::uint_len(ret->parsing_len, ret->trailing_char_bits), fp) != basics::uint_len(ret->parsing_len, ret->trailing_char_bits)){
        std::cerr<<"ERROR: load: trailing char "<<std::endl;
        delete ret; return NULL;
    }
    ret->phrases = DeltaCodes::load(fp);
    if(ret->phrases == NULL){
        std::cerr<<"ERROR: load: deltacodes "<<std::endl;
        delete ret; return NULL;
    }
    ret->wt_depth = static_sequence_wt::load(fp);
    if(ret->wt_depth == NULL){
        std::cerr<<"ERROR: load: depth WT"<<std::endl;
        return NULL; 
    }
    ret->perm = basics::static_permutation::load(fp);
    if(ret->perm == NULL){
        std::cerr<<"ERROR: load: permutation "<<std::endl;
        delete ret; return NULL;
    }
    ret->sources = DeltaCodes::load(fp);
    if(ret->sources == NULL){
        std::cerr<<"ERROR: load: deltacodes sources"<<std::endl;
        delete ret; return NULL;
    }

    return ret;
}

/*default constructor*/
static_compressed_text_lz77::static_compressed_text_lz77(){
        sigma_mapper = NULL;
        tlen = 0;
        parsing_len=0;
        trailing_char=NULL;
        trailing_char_bits=0;
        phrases=NULL;
        sources = NULL;
        perm = NULL;
        wt_depth = NULL;
        store_sst_ids = 1;
        binaryRev = 0;
        binarySst = 0;
}

/** Returns the substring of the text in a new allocated array*/
unsigned char* static_compressed_text_lz77::_display(unsigned int start, unsigned int end){
    unsigned char* answer = new unsigned char[end-start+2];//we add a trailing '\0' to the returned string
    answer[end-start+1] = '\0';
    _charextractLZ77(answer,start,end,0,end-start);
    return answer;
}

void static_compressed_text_lz77::_charextractLZ77(unsigned char* answer, unsigned start, unsigned end, unsigned int offset_start, unsigned int offset_end){
    while(1){//start<=end
        unsigned int last_phrase = 0;
        unsigned int r = this->phrases->rank_select(end,&last_phrase);
        if(last_phrase == end){
            answer[offset_end] = basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1;//the array doesnt store the 0
            if(start<end){
                end = end - 1;
                offset_end = offset_end - 1;
            }else{
                return;
            }
        }else{
            unsigned int source = _source(r);
            if(start>last_phrase){
                start = source + (start - last_phrase) -1;
                end = source + (end - last_phrase) -1;
                
            }else{
                _charextractLZ77(answer, source, source + (end - last_phrase) -1, offset_end - (end-last_phrase-1), offset_end);
                answer[offset_end - (end-last_phrase)] = basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1;//the array doesnt store the 0
                if(start<last_phrase){
                    offset_end = offset_end - (end-last_phrase) - 1;
                    end = last_phrase - 1;
                }else{
                    return;
                }
            }
        }
    }
}
size_t static_compressed_text_lz77::LCP( const std::string& P, size_t p_i, size_t t_j ){
    size_t n = length()-1;
    if( t_j >= n ) return 0;
    if( p_i >= P.size() ) return 0;
    size_t p_first_mismatch = m_find_mismatch_forward( P, p_i, t_j, n-1 );
    return p_first_mismatch - p_i;
}

size_t static_compressed_text_lz77::LCS( const std::string& P, size_t p_i, size_t t_j ){
    size_t n = length()-1;
    if( t_j >= n ) return 0;
    if( p_i >= P.size() ) return 0;
    std::pair<size_t,unsigned char> ret = m_find_matchspan_backward( P, p_i, 0, t_j );
    size_t p_last_match = ret.first;
    return p_i+1 - p_last_match;
}

std::pair<size_t, unsigned char> static_compressed_text_lz77::LCS_char( const std::string& P, size_t p_i, size_t t_j ){
    size_t n = length()-1;
    if( t_j >= n ) return std::make_pair(0,(unsigned char)-1);
    if( p_i >= P.size() ) return std::make_pair(0,(unsigned char)-1);
    std::pair<size_t,unsigned char> ret = m_find_matchspan_backward( P, p_i, 0, t_j );
    size_t p_last_match = ret.first;
    return std::make_pair( p_i+1 - p_last_match, ret.second );
}

unsigned char static_compressed_text_lz77::extract( unsigned i ) { 
    while(1) {
        unsigned int last_phrase = 0;
        unsigned int r = this->phrases->rank_select(i,&last_phrase);
        if(last_phrase == i){
            return sigma_mapper->unmapchar( basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1 ); //the array doesnt store the 0
        }else{
            unsigned int source = _source(r);
            i = source + (i - last_phrase) -1;
        }
    }
}

size_t static_compressed_text_lz77::m_find_mismatch_forward( const std::string& P, size_t p_i, size_t start, size_t end ){
    while( p_i < P.size() ){//start<=end
        unsigned int last_phrase = 0;
        unsigned int r = this->phrases->rank_select(start,&last_phrase);
        if(last_phrase == start){
            char c = sigma_mapper->unmapchar( basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1 ); //the array doesnt store the 0
            if( c != P[p_i] ) return p_i;
            p_i++;
            if(start<end){
                start ++;
            }else{
                return p_i;
            }
        }else{
            unsigned int source = _source(r);
            unsigned int next_phrase = this->phrases->select(r+1);
            if(end < next_phrase ){
                start = source + (start - last_phrase) -1;
                end = source + (end - last_phrase) -1;
                
            }else{
                size_t p_i_ = m_find_mismatch_forward( P, p_i, source + (start - last_phrase) - 1, source + (next_phrase - last_phrase) - 2 );
                if( p_i_ - p_i < next_phrase - start ) return p_i_;
                p_i = p_i_;

                char c = sigma_mapper->unmapchar( basics::get_field(this->trailing_char,this->trailing_char_bits,r)+1 );//the array doesnt store the 0
                if( c != P[p_i] ) return p_i;
                p_i++;
                if( end > next_phrase){
                    start = next_phrase+1;
                }else{
                    return p_i;
                }
            }
        }
    }
    return p_i;
}

std::pair<size_t,unsigned char> static_compressed_text_lz77::m_find_matchspan_backward( const std::string& P, size_t p_i, size_t start, size_t end ){
    while(1){//start<=end
        unsigned int last_phrase = 0;
        unsigned int r = this->phrases->rank_select(end,&last_phrase);
        if(last_phrase == end){
            char c = sigma_mapper->unmapchar( basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1 ); //the array doesnt store the 0
            if( P[p_i] != c ) return std::make_pair( p_i+1, (unsigned char)c );

            if(start>=end || p_i == 0 ) {
                return std::make_pair( p_i, (unsigned char)-1 );
            }
            end--;
            p_i--;
        }else{
            unsigned int source = _source(r);
            if(start>last_phrase){
                start = source + (start - last_phrase) -1;
                end = source + (end - last_phrase) -1;
                
            }else{
                std::pair<size_t, unsigned char> ret_ = m_find_matchspan_backward( P, p_i, source, source + (end - last_phrase) -1 );
                size_t p_i_ = ret_.first;
                if( p_i_ == 0 || p_i+1 - p_i_ < end - last_phrase ) return ret_;
                p_i = p_i_-1;

                char c = sigma_mapper->unmapchar( basics::get_field(this->trailing_char,this->trailing_char_bits,r-1)+1 );//the array doesnt store the 0
                if( P[p_i] != c ) return std::make_pair( p_i + 1, (unsigned char) c );

                if(start>=last_phrase || p_i == 0 ) {
                    return std::make_pair( p_i, (unsigned char)-1 );
                }
                end = last_phrase - 1;
                p_i --;
            }
        }
    }
}

int static_compressed_text_lz77::_source(unsigned int i){
    unsigned int r = perm->pi(i)+1;
    int pos = this->sources->select(r);
    return pos-r;
}
unsigned int static_compressed_text_lz77::_source_length(unsigned int i){
    //if(i==0)return 1;
    unsigned int length;
    phrases->select2(i,&length);
    return length;
}

}
