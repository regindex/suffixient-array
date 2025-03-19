#ifndef BITPACKED_TEXT_ORACLE_HPP
#define BITPACKED_TEXT_ORACLE_HPP

#include <common_.hpp>

namespace suffixient {

# define PACKED_INT_SIZE 2

class bitpacked_text_oracle
{
public:

    bitpacked_text_oracle(){};

    void build(std::string input_file_path)
    {
        std::ifstream file_text(input_file_path, std::ios::binary);
        file_text.seekg(0, std::ios::end);
        this->N = file_text.tellg();
        file_text.seekg(0, std::ios::beg);

        this->T.resize(this->N);
        for(usafe_t i=0;i<N;++i)
        {
            char_t c;
            file_text.read(reinterpret_cast<char*>(&c), sizeof(char_t));
            if(dna_to_code_table[c] > 3)
                { std::cerr << "Non DNA character detected!" << std::endl; exit(1); }
            T[i] = dna_to_code_table[c];
        }
    } 

    usafe_t store(std::string output_file_path)
    {
        usafe_t w_bytes = (this->N * PACKED_INT_SIZE)/8;

        return w_bytes;
    }

    void load(std::string input_file_path)
    {
        build(input_file_path);
    }

    usafe_t size(){ return (T.size() * PACKED_INT_SIZE)/8 + sizeof(N); }

    unsigned char extract(usafe_t i){ return code_to_dna_table[this->T[i]]; }
    
    usafe_t LCP(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((pattern.size()-p),(this->N-t));

        while(available_chars > 0)
        {
            if(pattern[p+matched_chars] != code_to_dna_table[this->T[t+matched_chars]])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }
    /*
    usafe_t LCP_(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((pattern.size()-p),(this->N-t));
        usafe_t available_blocks = available_chars / (64/PACKED_INT_SIZE);

        while(available_blocks > 0)
        {
            usafe_t T_int = this->T.get_int((t+matched_chars)*PACKED_INT_SIZE);
            usafe_t P_int = 0;
            set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&P_int), pattern, 
                             (64/PACKED_INT_SIZE), p+matched_chars, (64/PACKED_INT_SIZE));
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0){ return matched_chars + (__builtin_ctzll(xor_int) / PACKED_INT_SIZE); }

            matched_chars += 64/PACKED_INT_SIZE;
            available_chars -= 64/PACKED_INT_SIZE;
            available_blocks--;
        }
        {
            usafe_t T_int = this->T.get_int((t+matched_chars)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t P_int = 0;
            set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&P_int), pattern, available_chars, p+matched_chars, available_chars);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0){ return matched_chars + (__builtin_ctzll(xor_int) / PACKED_INT_SIZE); }

            matched_chars += available_chars;
        }

        return matched_chars;
    }
    usafe_t LCP(sdsl::int_vector<PACKED_INT_SIZE>& P, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((P.size()-p),(this->N-t));
        usafe_t available_blocks = available_chars / (64/PACKED_INT_SIZE);

        while(available_blocks > 0)
        {
            usafe_t T_int = this->T.get_int((t+matched_chars)*PACKED_INT_SIZE);
            usafe_t P_int = P.get_int((p+matched_chars)*PACKED_INT_SIZE);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0){ return matched_chars + (__builtin_ctzll(xor_int) / PACKED_INT_SIZE); }

            matched_chars += 64/PACKED_INT_SIZE;
            available_chars -= 64/PACKED_INT_SIZE;
            available_blocks--;
        }
        {
            usafe_t T_int = this->T.get_int((t+matched_chars)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t P_int = P.get_int((p+matched_chars)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0){ return matched_chars + (__builtin_ctzll(xor_int) / PACKED_INT_SIZE); }

            matched_chars += available_chars;
        }

        return matched_chars;
    }
    */
    usafe_t LCS(std::string& pattern, size_t p, size_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != code_to_dna_table[this->T[t-matched_chars]])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }
    /*
    std::pair<usafe_t,char_t> LCS_char_(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);
        usafe_t available_blocks = available_chars / (64/PACKED_INT_SIZE);

        while(available_blocks > 0)
        {
            usafe_t T_int = this->T.get_int((t-matched_chars-(64/PACKED_INT_SIZE)+1)*PACKED_INT_SIZE);
            usafe_t P_int = 0;
            set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&P_int), pattern, 
                             (64/PACKED_INT_SIZE), p-matched_chars-(64/PACKED_INT_SIZE)+1, (64/PACKED_INT_SIZE));
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0)
            {
                usafe_t lcs_ = matched_chars + (__builtin_clzll(xor_int) / PACKED_INT_SIZE);
                return std::make_pair(lcs_,code_to_dna_table[this->T[t-lcs_]]);
            }

            matched_chars += 64/PACKED_INT_SIZE;
            available_chars -= 64/PACKED_INT_SIZE;
            available_blocks--;
        }
        {
            usafe_t T_int = this->T.get_int((t-matched_chars)*PACKED_INT_SIZE-(available_chars-1)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t P_int = 0;
            set_uint_DNA_inv(reinterpret_cast<uint8_t*>(&P_int), pattern, available_chars, p-matched_chars-available_chars+1, available_chars);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0)
            {
                usafe_t lcs_ = matched_chars + ( (__builtin_clzll(xor_int) - ((64/PACKED_INT_SIZE)-available_chars)*PACKED_INT_SIZE) / PACKED_INT_SIZE);
                return std::make_pair(lcs_,code_to_dna_table[this->T[t-lcs_]]);
            }

            matched_chars += available_chars;
        }

        return std::make_pair(matched_chars,-1);
    }
    std::pair<usafe_t,char_t> LCS_char(sdsl::int_vector<PACKED_INT_SIZE>& P, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);
        usafe_t available_blocks = available_chars / (64/PACKED_INT_SIZE);

        while(available_blocks > 0)
        {
            usafe_t T_int = this->T.get_int((t-matched_chars-(64/PACKED_INT_SIZE)+1)*PACKED_INT_SIZE);
            usafe_t P_int = P.get_int((p-matched_chars-(64/PACKED_INT_SIZE)+1)*PACKED_INT_SIZE);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0)
            {
                usafe_t lcs_ = matched_chars + (__builtin_clzll(xor_int) / PACKED_INT_SIZE);
                return std::make_pair(lcs_,code_to_dna_table[this->T[t-lcs_]]);
            }

            matched_chars += 64/PACKED_INT_SIZE;
            available_chars -= 64/PACKED_INT_SIZE;
            available_blocks--;
        }
        {
            usafe_t T_int = this->T.get_int((t-matched_chars)*PACKED_INT_SIZE-(available_chars-1)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t P_int = P.get_int((p-matched_chars)*PACKED_INT_SIZE-(available_chars-1)*PACKED_INT_SIZE,available_chars*PACKED_INT_SIZE);
            usafe_t xor_int = T_int ^ P_int;

            if(xor_int != 0)
            {
                usafe_t lcs_ = matched_chars + ( (__builtin_clzll(xor_int) - ((64/PACKED_INT_SIZE)-available_chars)*PACKED_INT_SIZE) / PACKED_INT_SIZE);
                return std::make_pair(lcs_,code_to_dna_table[this->T[t-lcs_]]);
            }

            matched_chars += available_chars;
        }

        return std::make_pair(matched_chars,-1);
    }
    */
    std::pair<usafe_t,char_t> LCS_char(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != code_to_dna_table[this->T[t-matched_chars]])
                return std::make_pair(matched_chars,code_to_dna_table[this->T[t-matched_chars]]);

            matched_chars++;
            available_chars--;
        }

        return std::make_pair(matched_chars,-1);
    }
    
private:

   sdsl::int_vector<PACKED_INT_SIZE> T;  
   usafe_t N;
};

}  // namespace suffixient

#endif  // BITPACKED_TEXT_ORACLE_HPP