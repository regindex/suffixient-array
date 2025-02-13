#ifndef UNCOMPRESSED_TEXT_ORACLE_HPP
#define UNCOMPRESSED_TEXT_ORACLE_HPP

#include <common.hpp>

namespace suffixient {

class uncompressed_text_oracle
{
public:

    uncompressed_text_oracle(){};

    void build(std::string input_file_path)
    {
        std::ifstream file_text(input_file_path, std::ios::binary);
        file_text.seekg(0, std::ios::end);
        this->N = file_text.tellg();
        file_text.seekg(0, std::ios::beg);

        this->T.resize(this->N);
        file_text.read(reinterpret_cast<char*>(&(this->T[0])), this->N);
    }

    usafe_t store(std::string output_file_path)
    {
        usafe_t w_bytes = this->N;

        return w_bytes;
    }

    void load(std::string input_file_path)
    {
        build(input_file_path);
    }

    usafe_t size(){ return T.size() + sizeof(N); }

    unsigned char extract(usafe_t i){ return this->T[i]; }

    usafe_t LCP(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((pattern.size()-p),(this->N-t));

        while(available_chars > 0)
        {
            if(pattern[p+matched_chars] != this->T[t+matched_chars])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }

    usafe_t LCS(std::string& pattern, size_t p, size_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != this->T[t-matched_chars])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }

    std::pair<usafe_t,char_t> LCS_char(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);
        ////std::cout << "available chars: " << available_chars << std::endl;

        while(available_chars > 0)
        {
            ////std::cout << pattern[p-matched_chars] << " <***> "<< this->T[t-matched_chars] << std::endl;
            if(pattern[p-matched_chars] != this->T[t-matched_chars])
                return std::make_pair(matched_chars,this->T[t-matched_chars]);

            matched_chars++;
            available_chars--;
        }

        return std::make_pair(matched_chars,-1);
    }

private:

   sdsl::int_vector<8> T;  
   usafe_t N;
};

}  // namespace suffixient

#endif  // UNCOMPRESSED_TEXT_ORACLE_HPP
