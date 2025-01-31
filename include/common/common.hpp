// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 *  common.hpp: DESCRIPTION
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <sys/stat.h>

#ifndef M64
	#define M64 0
#endif
#if M64
    typedef uint64_t uint_t;
    typedef int64_t  int_t;
#else
    typedef uint32_t uint_t;
    typedef int32_t  int_t;
#endif
typedef int64_t safe_t;
typedef uint64_t usafe_t;
typedef bool bool_t;
typedef char char_t;
typedef uint8_t uchar_t;

#define DEF_BUFFER_SIZE 5000
#define SIGMA 128

static unsigned char dna_to_code_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 5
};

static unsigned char bit_mask_table[16] = {
	0,   0,   0,   0, 
	64,  16,  4,   1, 
	128, 32,  8,   2, 
	192, 48, 12,   3
}; 

inline bool_t contained(safe_t s1, safe_t e1, safe_t s2, safe_t e2)
{
    if( (s2 >= s1) and (e2 <= e1) )
        return true;
    return false;
}

void inline set_uint_DNA(uchar_t* t, std::string p, uchar_t size)
{
    size /= 8;
	for(uint_t i=0;i<p.size();++i)
		t[3-(i/4)] |= 
		bit_mask_table[(dna_to_code_table[p[i]]*4)+(i%4)];
}

usafe_t get_5bytes_uint(uchar_t *a, uint_t i=0)
{
  uint_t offset = (i+1)*5-1;
  usafe_t ai = 0;
  for(uint_t j=0;j<5;j++) ai = (ai << 8) | a[offset-j];
    
  return ai;
}

template<typename T>
void read_file(const char *filename, std::vector<T>& ptr)
{
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        std::cerr << "open() file " + std::string(filename) + " failed\n";

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        std::cerr << "stat() file " + std::string(filename) + " failed\n";

    if(filestat.st_size % sizeof(T) != 0)
        std::cerr << "invilid file " + std::string(filename) + "\n";

    size_t length = filestat.st_size / sizeof(T);
    ptr.resize(length);

    if ((fread(&ptr[0], sizeof(T), length, fd)) != length) 
        std::cerr << "fread() file " + std::string(filename) + " failed\n";

    fclose(fd);
}

#endif // COMMON_HPP_