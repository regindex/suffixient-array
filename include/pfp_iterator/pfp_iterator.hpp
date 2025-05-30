// Copyright (c) 2024, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.
/*
    Adapted from a the file "pfp_lcp.hpp" in "https://github.com/maxrossi91/pfp-thresholds"
    by Missimiliano Rossi 2020.
*/

#ifndef _PFP_ITERATOR_HH
#define _PFP_ITERATOR_HH

//#include <common.hpp>
#include "common.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C"
{
#include <gsacak.h>
}

//#include <pfp.hpp>
//#include <priority_queue.hpp>
#include "pfp.hpp"
#include "priority_queue.hpp"

class pfp_iterator{
public:

    pfp_iterator(pf_parsing &pfp_, std::string filename) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                head(0)
    { assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict); } 

    void process_next_block()
    {
        while (not is_finished())
        {
            if(is_valid(curr)){
                // Compute the next character of the BWT of T
                // Store the list of all phrase ids with the same suffix.
                std::vector<phrase_suffix_t> same_suffix(1, curr); 

                next = curr;

                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }
                }

                // compute lcp suffix length
                lcp_suffix = compute_lcp_suffix(curr,prev);

                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                }

                assert(pq.size() > 0);

                // update prev
                prev = same_suffix.back();
                first = true;
                
                break;
            }
            else
            {
                inc(curr);
            }
        }
    }

    //void next_stream_pos()
    bool operator ++()
    {
        if(pq.empty())
        {
            curr = next;
            process_next_block();

            if( pq.empty() )
                return false;
        }

        auto curr_occ = pq.top();
        pq.pop();

        if (!first)
        {
            // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
            lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
        }
        first = false;

        // update LCP, SA and BWT entries
        update_lcp(lcp_suffix);

        update_sa(curr, *curr_occ.first);

        update_bwt(curr_occ.second.second);

        // Update prevs
        prev_occ = *curr_occ.first;

        // Update pq
        curr_occ.first++;
        if (curr_occ.first != curr_occ.second.first)
            pq.push(curr_occ);

        j += 1;

        return true;
    }

    bool is_finished()
    {
        // if the queue is empty and there are no more BWT blocks we stop
        return ( pq.empty() and (curr.i >= pf.dict.saD.size()) );
    }

    size_t get_sa()
    {
        return sas;
    }

    uint8_t get_bwt()
    {
        return head;
    }

    int_t get_lcp()
    {
        return lcpe;
    }

private:

    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        int_da sn = 0;
        uint8_t bwt_char = 0;
    } phrase_suffix_t;

    pf_parsing& pf;
    std::vector<size_t> min_s;  // Value of the minimum lcp_T in each run of BWT_T
    std::vector<size_t> pos_s;  // Position of the minimum lcp_T in each run of BWT_T

    phrase_suffix_t next;
    phrase_suffix_t curr;
    phrase_suffix_t prev;

    uint8_t head;
    size_t length = 0; // Length of the current run of BWT_T
    size_t j = 0;
    size_t sas = 0; // suffix array sample
    int_t lcpe = 0;

    size_t prev_occ;
    bool first;
    int_t lcp_suffix;

    // define the priority queue data structure
    typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;
    PriorityQueue<pq_t> pq;

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn);
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.pars.ilist.size()));
        s.suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        if(is_valid(s))
            s.bwt_char = (s.sn == pf.w ? 0 : pf.dict.d[s.sn - 1]);
        return true;
    }

    inline bool is_valid(phrase_suffix_t& s)
    {
        // avoid the extra w # at the beginning of the text
        if (s.sn < pf.w)
            return false;
        // Check if the suffix has length at least w and is not the complete phrase.
        if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
            return false;
        
        return true;
    }
    
    inline int_t min_s_lcp_T(size_t left, size_t right)
    {
        // assume left < right
        if (left > right)
            std::swap(left, right);

        assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

        return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
    }

    inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev)
    {
        int_t lcp_suffix = 0;

        if (j > 0)
        {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[curr.i];
            for (size_t k = prev.i + 1; k < curr.i; ++k)
            {
                lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= curr.suffix_length && curr.suffix_length == prev.suffix_length)
            {
                // Compute the minimum s_lcpP of the phrases following the two phrases
                // we take the first occurrence of the phrase in BWT_P
                size_t left = pf.pars.ilist[pf.pars.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.pars.ilist[pf.pars.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
                lcp_suffix += min_s_lcp_T(left,right);
            }
        }

        return lcp_suffix;
    }

    inline void update_sa(phrase_suffix_t &curr, size_t pos)
    {
        sas = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(sas < (pf.n - pf.w + 1ULL));
    }

    inline void update_bwt(uint8_t next_char)
    {
        if (head != next_char)
            head = next_char;
    }

    inline void update_lcp(int_t val)
    {
        lcpe = val;
    }

};

#endif /* end of include guard: _PFP_ITERATOR_HH */
