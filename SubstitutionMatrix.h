//
//  SubstitutionMatrix.h
//  AlgoUmMilotZuZerstoeren
//
//  Created by Martin Steinegger on 23.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef AlgoUmMilotZuZerstoeren_SubstitutionMatrix_h
#define AlgoUmMilotZuZerstoeren_SubstitutionMatrix_h
#include <map>

class SubstitutionMatrix
{
public:    
    SubstitutionMatrix(float ** matrix,std::map<char,int> rowcolLookup);
    ~SubstitutionMatrix();
    char *short2aa;
    short *aa2short;
    
    unsigned char * get_aa_substitution_vector(short aa);
    unsigned short get_substituion_score(short aa_q,short aa_t);
    unsigned short get_substituion_score(char aa_q,char aa_t);
    void print();
    short scal_factor;
    short middle_of_rang;
private:
    size_t size_rowcol;
    unsigned char ** substitution_matrix;
    
};


#endif
