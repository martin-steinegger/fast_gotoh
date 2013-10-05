//
//  SequenceFour.h
//  AlgoUmMilotZuZerstoeren
//
//  Created by Martin Steinegger on 22.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef AlgoUmMilotZuZerstoeren_Sequence8_h
#define AlgoUmMilotZuZerstoeren_Sequence8_h

#include <tmmintrin.h>

class SequenceEight
{
public:
    SequenceEight(size_t maxres,short* aa2short,char* short2aa);
    ~SequenceEight();
    size_t L;                    // length of sequence
    __m128i * sequence; 
    size_t sequence_count;
    short *sequence_short;
    short *aa2short;
    char  *short2aa ;
    size_t * seq_size_array;
    // Map
    void MapFourSEQsToSEQ8(const char * seq1,const char * seq2,const char * seq3,const char * seq4,
                           const char * seq5,const char * seq6,const char * seq7,const char * seq8);
    void MapOneSEQToSEQ8(const char *seq);
    void MapSEQArrayToSEQ8(const char ** seq,size_t element_count);
    // Access
    short get_first_aa_short_index(size_t pos);
    size_t get_sequence_count();
    size_t get_sequence_size(size_t pos);
    char get_char_by_pos(size_t pos,size_t elem_index);
    short get_aa_short_index(size_t pos,size_t elem_index);
    // Debug
    void print();
private: 
    const char ** sequence_array;
    
};


#endif
