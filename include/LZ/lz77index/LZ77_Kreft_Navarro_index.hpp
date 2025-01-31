#ifndef LZ77_DISPLAY_HPP
#define LZ77_DISPLAY_HPP

#include <static_selfindex.h>
#include <utils.h>

namespace lz77 {

//template <typename Value, const Value EMPTY_VALUE = -1>
class LZ77_Kreft_Navarro_index
{
public:

    LZ77_Kreft_Navarro_index(){};

    void build(std::string input_file, size_t window_size_ = 8)
    {
        this->window_size = window_size_;
        std::string index_file = input_file + ".lz77";
        idx = lz77index::static_selfindex_lz77::build(&input_file[0],&index_file[0], br, bs, ss);
    }

    size_t LCP(std::string& pattern, size_t p, size_t t)
    {
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
                    return matched_chars + i;
            }

            matched_chars += i;
            p += i; t += i;
            available_chars -= i;
            //window_size_ *= window_size;
            window_size_ *= 2;
        }

        return matched_chars;
    }

    size_t LCS(std::string& pattern, size_t p, size_t t)
    {
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
                    return matched_chars + i;
            }

            matched_chars += i;
            p -= i; t -= i;
            available_chars -= i;
            //window_size_ *= window_size;
            window_size_ *= 2;
        }

        return matched_chars;
    }

private:

   lz77index::static_selfindex* idx;  
   size_t window_size;
   unsigned char br = 0, bs = 0, ss = 0;
};

}  // namespace lz77

#endif  // LZ77_DISPLAY_HPP
