//
//
//  Created by Martin Steinegger
//  Copyright (c) 2012 -. All rights reserved.
//
#include "SequenceEight.h"
#include "util.h"
#include <algorithm>
#include <float.h>    // FLT_MIN
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <string.h>


/////////////////////////////////////////////////////////////////////////////////////
// Object constructor
/////////////////////////////////////////////////////////////////////////////////////
const int VECTOR_SIZE=8;

SequenceEight::SequenceEight(size_t maxres,short* aa2short,char* short2aa)
{
    this->sequence = (__m128i *)memalign(16,VECTOR_SIZE*maxres*sizeof(short)); 
    this->sequence_short = (short*)memalign(16,VECTOR_SIZE*maxres*sizeof(short)); 
    memset (sequence_short,(unsigned short)0,VECTOR_SIZE*maxres);
    memset (sequence,(unsigned short)0,VECTOR_SIZE*maxres);
    
    this->aa2short = aa2short;
    this->short2aa = short2aa;
    this->seq_size_array= new size_t[VECTOR_SIZE]; 
    
}



/////////////////////////////////////////////////////////////////////////////////////
// Object destructor
/////////////////////////////////////////////////////////////////////////////////////
SequenceEight::~SequenceEight()
{
    free(sequence); 
}

void SequenceEight::MapOneSEQToSEQ8(const char *seq){
    MapFourSEQsToSEQ8(seq,seq,seq,seq,seq,seq,seq,seq);
}

void SequenceEight::MapFourSEQsToSEQ8(const char * seq1,const char * seq2,const char * seq3,const char * seq4,
                                      const char * seq5,const char * seq6,const char * seq7,const char * seq8){
    const char * seqarr[VECTOR_SIZE]={seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8};  
    MapSEQArrayToSEQ8(seqarr, VECTOR_SIZE);
}

void SequenceEight::MapSEQArrayToSEQ8(const char ** seq,size_t element_count){
    this->L=0;
    this->sequence_count=element_count;
    for(size_t i = 0 ; i<element_count; i++){
        size_t curr_size = strlen(seq[i]);
        this->L=std::max(this->L,curr_size);     
        this->seq_size_array[i]=curr_size;
    }
    short * sequence_scalar=(short *)sequence;
    
    for(size_t seq_i = 0; seq_i < element_count; seq_i++){
        const char * curr = seq[seq_i];
        size_t curr_len=strlen(curr);
        for(size_t i = 0; i < curr_len; i++){
            short aa_short=aa2short[curr[i]];
            sequence_scalar[i*VECTOR_SIZE+seq_i] = aa_short;
            sequence_short [i*VECTOR_SIZE+seq_i] = aa_short;
        }
        for(size_t i = curr_len; i < this->L; i++){
            sequence_scalar[i*VECTOR_SIZE+seq_i] = 0;
            sequence_short [i*VECTOR_SIZE+seq_i] = 0;
            
        }
    }
    
}

short SequenceEight::get_first_aa_short_index(size_t pos){
    return sequence_short[pos*VECTOR_SIZE];
}


short SequenceEight::get_aa_short_index(size_t pos,size_t elem_index){
    return sequence_short[pos*VECTOR_SIZE+elem_index];
}

size_t SequenceEight::get_sequence_count(){
    return this->sequence_count;
}

size_t SequenceEight::get_sequence_size(size_t pos){
    return seq_size_array[pos];
}

char SequenceEight::get_char_by_pos(size_t pos,size_t elem_index){
    return short2aa[get_aa_short_index(pos,elem_index)];
}


void SequenceEight::print() {
	for (int i = 0; i < this->L; i++) {
        std::cout << " ( ";
        for(int j = 0; j < VECTOR_SIZE; j++){
            short toprint= _mm_extract_epi16(this->sequence[i],j);
			std::cout << toprint << ",";
        }
        std::cout << " ) ";
	}
	std::cout << std::endl;	 
}

