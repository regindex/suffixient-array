#ifndef LZ77_COMPRESSED_TEXT_HPP
#define LZ77_COMPRESSED_TEXT_HPP

#include <static_compressed_text.h>
#include <utils.h>

namespace lz77 {

//template <typename Value, const Value EMPTY_VALUE = -1>
class LZ77_compressed_text
{
public:

    LZ77_compressed_text(){};

    void build(std::string input_file, size_t window_size_ = 8)
    {
        this->window_size = window_size_;
        std::string index_file = input_file + ".lz77";
        idx = lz77index::static_compressed_text_lz77::build(&input_file[0],&index_file[0], br, bs, ss);
    }

    void load(std::string index_file)
    {
        FILE* lz77index_fp = fopen(index_file.c_str(),"r");
        this->idx = lz77index::static_compressed_text_lz77::load(lz77index_fp);
        fclose(lz77index_fp);
        this->window_size = 8;
    }

    size_t size(){ return idx->size(); }

    unsigned char extract(size_t i)
    {
        return idx->extract(i);
        /*      
        unsigned char* window = idx->display(i, i+1);

        unsigned char c = *window;
        delete[] window;
        
        return c;
        */
    }

    size_t LCP(std::string& pattern, size_t p, size_t t)
    {
        return idx->LCP( pattern, p, t );
        /*
        size_t matched_chars = 0;
        size_t window_size_ = this->window_size;
        size_t available_chars = std::min((pattern.size()-p),(idx->length()-t));

        while(available_chars > 0)
        {
            window_size_ = std::min(window_size_,available_chars);
            unsigned char* window = idx->display(t, t+window_size_);

            size_t i=0;
            for(;i<window_size_;++i)
            {
                if(*(window+i) != pattern[p+i])
                {
                    delete[] window;
                    
                    return matched_chars + i;
                }
            }

            matched_chars += i;
            p += i; t += i;
            available_chars -= i;
            //window_size_ *= window_size;
            window_size_ *= 2;

            delete[] window;
        }

        return matched_chars;
        */
    }

    size_t LCS(std::string& pattern, size_t p, size_t t)
    {
        return idx->LCS( pattern, p, t );
        /*
        size_t matched_chars = 0;
        size_t window_size_ = this->window_size;
        size_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            window_size_ = std::min(window_size_,available_chars);
            unsigned char* window = idx->display(t-window_size_+1, t+1);

            size_t i=0;
            for(;i<window_size_;++i)
            {
                if(*(window+(window_size_-1-i)) != pattern[p-i])
                {
                    delete[] window;

                    return matched_chars + i;
                }
            }

            matched_chars += i;
            p -= i; t -= i;
            available_chars -= i;
            //window_size_ *= window_size;
            window_size_ *= 2;

            delete[] window;
        }

        return matched_chars;
        */
    }

    std::pair<size_t,unsigned char> LCS_char(std::string& pattern, size_t p, size_t t,
                                    bool adaptable_window = false)
    {
        return idx->LCS_char( pattern, p, t );
        /*
        size_t matched_chars = 0;
        size_t available_chars = std::min(p+1,t+1);
        size_t window_size_ = std::min(this->window_size,available_chars);
        // size_t window_size_ = available_chars;

        while(available_chars > 0)
        {
            window_size_ = std::min(window_size_,available_chars);
            unsigned char* window = idx->display(t-window_size_+1, t+1);

            size_t i=0;
            for(;i<window_size_;++i)
            {
                if(*(window+(window_size_-1-i)) != pattern[p-i])
                {
                    unsigned char c = *(window+(window_size_-1-i));
                    delete[] window;

                    return std::make_pair(matched_chars + i, c);
                }
            }

            matched_chars += i;
            p -= i; t -= i;
            available_chars -= i;
            window_size_ *= 2;

            delete[] window;
        }

        return std::make_pair(matched_chars,-1);
        */
    }

private:

   lz77index::static_compressed_text_lz77* idx;  
   size_t window_size;
   unsigned char br = 0, bs = 0, ss = 0;
};

}  // namespace lz77

#endif  // LZ77_COMPRESSED_TEXT_HPP
