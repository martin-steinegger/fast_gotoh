//
//  GotohEight.h
//  AlgoUmMilotZuZerstoeren
//
//  Created by Martin Steinegger on 23.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#ifndef AlgoUmMilotZuZerstoeren_GotohEight_h
#define AlgoUmMilotZuZerstoeren_GotohEight_h
#include "SequenceEight.h"
#include "SubstitutionMatrix.h"
#include <tmmintrin.h>
#include <vector>

const int VEC_SIZE=8;
const int ALI_TYPE_GLOBAL =1;
const int ALI_TYPE_LOCAL =2;
const int ALI_TYPE_FREESHIFT =3;

class GotohEight
{
public:
    struct GotohMatrix {
        short ** scoreIns;
        short ** scoreDel;
        short ** scoreMat;
        size_t i[VEC_SIZE];
        size_t j[VEC_SIZE];
        float score[VEC_SIZE];
    };
    
    GotohEight(size_t maxres,SubstitutionMatrix * subMatrix,short go,short ge,int type);
    ~GotohEight();
    GotohMatrix calc_matrix(SequenceEight* q, SequenceEight* t);
    std::vector<char *> backtrace(GotohEight::GotohMatrix matrix,SequenceEight* q, SequenceEight* t);
    float calc_check_score(char *q, char * t);
private:
    size_t maxres;
    int ali_type;
    SubstitutionMatrix * subMatrix;
    void init_global();
    void init_local();
    void init_freeshift();
    GotohMatrix global(SequenceEight* q, SequenceEight* t);
    GotohMatrix freeshift(SequenceEight* q, SequenceEight* t);
    GotohMatrix local(SequenceEight* q, SequenceEight* t);
    __m128i create_profile(__m128i * score_matrix_vec01,__m128i * score_matrix_vec16,__m128i * template_sequence);
    
    GotohMatrix returnMatrix;
    short ** scoreIns;
    short ** scoreDel;
    short ** scoreMat;
    short go;
    short ge;
    short goe;
    short NEGATIV_INF;
    __m128i sixteen_vec;
    __m128i fiveteen_vec;
    char ** t_backtrace_seq_str;
    char ** q_backtrace_seq_str;
};
void print_m128_vec(__m128 toprint_vec);
#endif
