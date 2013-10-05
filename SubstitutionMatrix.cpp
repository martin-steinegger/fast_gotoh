//
//  SubstitutionMatrix.cpp
//  AlgoUmMilotZuZerstoeren
//
//  Created by Martin Steinegger on 23.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include "SubstitutionMatrix.h"
#include <iostream>
#include "util.h"
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string.h>
SubstitutionMatrix::SubstitutionMatrix(float ** matrix,
                                       std::map<char,int> rowcolLookup){
    
    this->size_rowcol=rowcolLookup.size();
    std::map<char, int>::iterator it = rowcolLookup.begin();
    int maxCharIndex=0;
    short2aa = new char[size_rowcol];
    for( ; it != rowcolLookup.end(); ++it )
    {
        char aa_char = it->first;
        int matrix_index = it->second+1;
        if(matrix_index > 15)
            matrix_index +=1;
        short2aa[matrix_index] = (short) aa_char;
        maxCharIndex=std::max(maxCharIndex,(int)aa_char);
    }       
    
    aa2short = new short[maxCharIndex+1];
    // Null should be a AminoAcid value
    for(size_t i=1; i<=size_rowcol; ++i){
        size_t index_of_aa =i;
        if(i > 15)
            index_of_aa += 1;
        aa2short[short2aa[index_of_aa]] = (short) index_of_aa;
    }
    // +2 for the 2 zero lines
    substitution_matrix=(unsigned char**)memalign(16,(this->size_rowcol+2)*sizeof(unsigned char *)); 
    substitution_matrix[0] =(unsigned char*)memalign(16,32*sizeof(unsigned char));
    substitution_matrix[16] =(unsigned char*)memalign(16,32*sizeof(unsigned char));
    memset (substitution_matrix[0],(unsigned char)0,32);
    memset (substitution_matrix[16],(unsigned char)0,32);

    
    float matrix_max_float=0.0f;
    bool isFloat = false;
    scal_factor=1;
    for(int row = 0; row < size_rowcol; row++){
        for(int col =0; col < size_rowcol; col++){
            matrix_max_float=std::max(matrix_max_float, matrix[row][col]);
            if(ceilf(matrix[row][col]) != matrix[row][col]){
                isFloat=true;
                scal_factor=10;
            }
        }
    }
    middle_of_rang=255.0f-(matrix_max_float)*((float)scal_factor);
    int row_minus=1;
    for(int row=1; row < size_rowcol+2; row++){
        substitution_matrix[row] =(unsigned char*)memalign(16,32*sizeof(unsigned char));
        memset (substitution_matrix[row],(unsigned char)0,32);
        
        if(row == 16){
            row_minus = 2;
            continue;
        }
        
        for(int i=0; i<2; i++){
            size_t to_index=(size_t)(16*(i+1));
            size_t from_index= (16*i);
            size_t matrix_residue = std::min(to_index, size_rowcol);
            for(size_t col=from_index; col <= matrix_residue+i; col++){
                if(col==from_index){
                    substitution_matrix[row][col]=(unsigned char)0;
                    continue;
                }
                int col_minus=(1*(i+1));
                substitution_matrix[row][col]=(unsigned char)(middle_of_rang + ((float)scal_factor)*matrix[row-row_minus][col-col_minus]);
            }
 
        }    
        
    }
}

SubstitutionMatrix::~SubstitutionMatrix(){
    free(short2aa);
    free(aa2short);
}

unsigned short SubstitutionMatrix::get_substituion_score(short aa_q,short aa_t){
    return substitution_matrix[aa_q][aa_t];
}

unsigned short SubstitutionMatrix::get_substituion_score(char aa_q,char aa_t){
    return substitution_matrix[aa2short[aa_q]][aa2short[aa_t]];
}


unsigned char * SubstitutionMatrix::get_aa_substitution_vector(short aa){
    
    return substitution_matrix[aa];
}


void SubstitutionMatrix::print() {
    for (int i = 0; i < size_rowcol+2; i++) {
        printf("\t%c",(char)short2aa[i]);
    }
    std::cout << std::endl;

	for (int i = 0; i < size_rowcol+2; i++) {
        printf("%c\t",(char)short2aa[i]);

		for (int j = 0; j < 32; j++) {
            printf("%u\t", substitution_matrix[i][j]);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;	 
}
