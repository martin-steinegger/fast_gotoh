#include "GotohEight.h"
#include <limits.h>
#include <iostream>
#include "util.h"
#include <stdio.h>
#include <string.h>



GotohEight::GotohEight(size_t maxres,SubstitutionMatrix * subMatrix,short go,short ge,int type){
    this->maxres = maxres;
    this->subMatrix = subMatrix;
    this->ali_type = type;
    this->go = go*subMatrix->scal_factor;
    this->ge = ge*subMatrix->scal_factor;
    this->goe = this->go+this->ge;
    this->NEGATIV_INF=SHRT_MIN-goe;
    //    this->backTraceMatch=(short**)memalign(16,(maxres+1)*sizeof(short *)); 
    this->scoreIns=(short**)memalign(16,((maxres+1))*sizeof(short *)); 
    this->scoreDel=(short**)memalign(16,((maxres+1))*sizeof(short *)); 
    this->scoreMat=(short**)memalign(16,((maxres+1))*sizeof(short *)); 
    for(int row=0; row <= maxres; row++){
        this->scoreIns[row] =(short*)memalign(16,(maxres+1)*VEC_SIZE*sizeof(short));
        this->scoreDel[row] =(short*)memalign(16,(maxres+1)*VEC_SIZE*sizeof(short));
        this->scoreMat[row] =(short*)memalign(16,(maxres+1)*VEC_SIZE*sizeof(short));
    }
    for(int row=0; row <= maxres; row++){
        for(int col=0; col <= maxres; col++){
            this->scoreIns[row][col] =0;
            this->scoreDel[row][col] =0;
            this->scoreMat[row][col] =0;
        }
    }
    switch(ali_type) {
        case ALI_TYPE_LOCAL:
            init_local();
            break;
        case ALI_TYPE_GLOBAL:
            init_global();
            break;
        case ALI_TYPE_FREESHIFT:
        default:
            init_local();
            break;
    }    

   

    
    
    this->sixteen_vec = _mm_set1_epi8(16);
    this->fiveteen_vec = _mm_set1_epi8(15);
    
    this->t_backtrace_seq_str = (char**)memalign(16,VEC_SIZE*sizeof(char *)); 
    this->q_backtrace_seq_str = (char**)memalign(16,VEC_SIZE*sizeof(char *)); 
    for(size_t i = 0; i < VEC_SIZE;i++){
        this->t_backtrace_seq_str[i]=(char*)memalign(16,(2*maxres+2)*sizeof(char));
        this->q_backtrace_seq_str[i]=(char*)memalign(16,(2*maxres+2)*sizeof(char));
    }
}


void GotohEight::init_global(){
    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    // Initialization of top row, i.e. cells (0,j)
    for (int j=0; j<=maxres; ++j) 
    {
        short cost=0;
        if(j!=0)
            cost=go+ge*(j);
        sMat_vec[0][j] = _mm_set1_epi16(cost);
        sIns_vec[0][j] = _mm_set1_epi16(NEGATIV_INF); 
    }
    // Initialize cells        
    // initialize at (i,0)
    for (int i=1; i<=maxres; ++i){ 
        short cost=go+ge*(i);
        sMat_vec[i][0]  = _mm_set1_epi16(cost);
        sIns_vec[i][0] = _mm_set1_epi16(cost);
        sDel_vec[i][0] = _mm_set1_epi16(NEGATIV_INF); 
    }
    
}



void GotohEight::init_freeshift(){
    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    // Initialization of top row, i.e. cells (0,j)
    for (int j=0; j<=maxres; ++j) 
    {
        
        sMat_vec[0][j] = _mm_set1_epi16(0);
        sIns_vec[0][j] = _mm_set1_epi16(NEGATIV_INF); 
        sDel_vec[0][j] = _mm_set1_epi16(NEGATIV_INF); 

    }
    // Initialize cells        
    // initialize at (i,0)
    for (int i=0; i<=maxres; ++i){ 
        sMat_vec[i][0]  = _mm_set1_epi16(0);
        sIns_vec[i][0] = _mm_set1_epi16(NEGATIV_INF);
        sDel_vec[i][0] = _mm_set1_epi16(NEGATIV_INF); 
    }
    
}


void GotohEight::init_local(){
    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    // Initialization of top row, i.e. cells (0,j)
    for (int j=0; j<=maxres; ++j) 
    {
        sMat_vec[0][j] = _mm_set1_epi16(0);
        sIns_vec[0][j] = _mm_set1_epi16(NEGATIV_INF);
        sDel_vec[0][j] = _mm_set1_epi16(NEGATIV_INF); 

    }
    // Initialize cells        
    // initialize at (i,0)
    for (int i=0; i<=maxres; ++i){ 
        sMat_vec[i][0]  = _mm_set1_epi16(0);
        sIns_vec[i][0] = _mm_set1_epi16(NEGATIV_INF);
        sDel_vec[i][0] = _mm_set1_epi16(NEGATIV_INF); 
    }
    
}



GotohEight::~GotohEight(){
    for(int row=0; row <= maxres; row++){
        free(scoreIns[row]);
        free(scoreMat[row]);
        free(scoreDel[row]);
    }
    free(scoreIns);
    free(scoreMat);
    free(scoreDel);
}



inline __m128i GotohEight::create_profile(__m128i * score_matrix_vec01
                                          ,__m128i * score_matrix_vec16
                                          ,__m128i * template_sequence){
    
    
    __m128i template01=_mm_packs_epi16(*template_sequence,*template_sequence);    
    // create slice mask
    // Example:
    //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
    //                      if lt 16 
    //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
    __m128i lookup_mask01=_mm_cmplt_epi8(template01,sixteen_vec);
    __m128i lookup_mask16=_mm_cmpgt_epi8(template01,fiveteen_vec);
    // slice index
    // Example:
    //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
    //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
    //                          min	
    //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
    __m128i lookup_index01=_mm_min_epu8(lookup_mask01,template01);
    __m128i lookup_index16=_mm_min_epu8(lookup_mask16,template01);
    // 2xmal array lookup
    __m128i score01=_mm_shuffle_epi8(*score_matrix_vec01,lookup_index01);
    __m128i score16 =_mm_shuffle_epi8(*score_matrix_vec16,lookup_index16);
    
    __m128i score_vec_8bit=_mm_add_epi8(score01,score16);
    __m128i score_vec_16bit = _mm_srli_epi16(_mm_unpacklo_epi8(score_vec_8bit, score_vec_8bit), 8);
    
    return score_vec_16bit;
}


GotohEight::GotohMatrix GotohEight::freeshift(SequenceEight* q, SequenceEight* t)
{
    

    // Variable declarations
    const __m128i goe_vec= _mm_set1_epi16(goe);
    const __m128i ge_vec = _mm_set1_epi16(ge);
    const __m128i score_middle_vec=_mm_set1_epi16(subMatrix->middle_of_rang);
    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    
    int i,j;      //query and template match state indices
    
    for (i=1; i<=q->L; ++i) // Loop through query positions i
    {
        
        // get score for query aa
        //        std::cout <<(signed short) _mm_extract_epi16(sMat_vec[i][0],0) << '\t';
        short query_aa_index=q->get_first_aa_short_index(i-1);
        __m128i * sub_score_vec=(__m128i *)subMatrix->get_aa_substitution_vector(query_aa_index);
        
        __m128i score_matrix_vec01=*sub_score_vec;
        __m128i score_matrix_vec16=*(sub_score_vec+1);  
        size_t t_size=t->L;
        for (j=1; j<=t_size; ++j) // Loop through template positions j
        {
            
            // calculate amino acid profile-profile scores
            // __m128i score_vec=create_profile(&score_matrix_vec01,&score_matrix_vec16,&t->sequence[j-1]);
            __m128i template01=_mm_packs_epi16(t->sequence[j-1],t->sequence[j-1]);    
            // create slice mask
            // Example:
            //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
            //                      if lt 16 
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            __m128i lookup_mask01=_mm_cmplt_epi8(template01,sixteen_vec);
            __m128i lookup_mask16=_mm_cmpgt_epi8(template01,fiveteen_vec);
            // slice index
            // Example:
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11  
            //                          min	
            //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
            __m128i lookup_index01=_mm_min_epu8(lookup_mask01,template01);
            __m128i lookup_index16=_mm_min_epu8(lookup_mask16,template01);
            // 2xmal array lookup
            __m128i score01=_mm_shuffle_epi8(score_matrix_vec01,lookup_index01);
            __m128i score16 =_mm_shuffle_epi8(score_matrix_vec16,lookup_index16);
            
            __m128i score_vec_8bit=_mm_add_epi8(score01,score16);
            __m128i score_vec_16bit = _mm_srli_epi16(_mm_unpacklo_epi8(score_vec_8bit, score_vec_8bit), 8);
            
            
            
            
            score_vec_16bit = _mm_sub_epi16(score_vec_16bit,score_middle_vec);
            sDel_vec[i][j] = _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i][j-1],goe_vec),
                                           _mm_add_epi16(sDel_vec[i][j-1],ge_vec)
                                           );
            sIns_vec[i][j] = _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i-1][j],goe_vec),
                                           _mm_add_epi16(sIns_vec[i-1][j],ge_vec)
                                           );
            
            
            sMat_vec[i][j] = _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i-1][j-1],score_vec_16bit), 
                                           _mm_max_epi16(sDel_vec[i][j],sIns_vec[i][j])
                                           );
            
//                       std::cout << (signed short) _mm_extract_epi16(sMat_vec[i][j],0) << '\t';
            
            
        } //end for j
//                std::cout << std::endl;	         
        
    } // end for i
    
    returnMatrix.scoreDel = this->scoreDel;
    returnMatrix.scoreMat = this->scoreMat;
    returnMatrix.scoreIns = this->scoreIns;
    for(size_t seq_index=0; seq_index<t->get_sequence_count(); seq_index++){
        size_t i=q->get_sequence_size(seq_index);
        size_t j=t->get_sequence_size(seq_index);
        short max = NEGATIV_INF;
        for(size_t row_index = i;row_index != 0;row_index--){
            short mat_col_score=scoreMat[row_index][j*(VEC_SIZE)+seq_index];
            if(max < mat_col_score){
                max = mat_col_score;
                returnMatrix.i[seq_index] = row_index;
                returnMatrix.j[seq_index] = j;

            }
        }
        for(size_t col_index = j;col_index != 0;col_index--){
            short mat_col_score=scoreMat[i][col_index*(VEC_SIZE)+seq_index];
            if(max < mat_col_score){
                max = mat_col_score;
                returnMatrix.j[seq_index] = col_index;
                returnMatrix.i[seq_index] = i;
            }
        }
        
        returnMatrix.score[seq_index] =(float) (max /(float)subMatrix->scal_factor);
    }
    
    
    
    return returnMatrix;
}
////////////////////////////////////////////////////////////////////////




GotohEight::GotohMatrix GotohEight::local(SequenceEight* q, SequenceEight* t)
{
    const __m128i goe_vec= _mm_set1_epi16(goe);
    const __m128i ge_vec = _mm_set1_epi16(ge);
    const __m128i score_middle_vec=_mm_set1_epi16(subMatrix->middle_of_rang);
    const __m128i null_vec=_mm_set1_epi16(0);
    const __m128i minus_one_vec =_mm_set1_epi16(-1);

    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    __m128i max_scores=_mm_set1_epi16(NEGATIV_INF);
    __m128i col_max_score=_mm_set1_epi16(NEGATIV_INF);
    __m128i i_pos_max=_mm_set1_epi16(0);
    int i,j;      //query and template match state indices
    
    for (i=1; i<=q->L; ++i) // Loop through query positions i
    {
        
        // get score for query aa
        //        std::cout <<(signed short) _mm_extract_epi16(sMat_vec[i][0],0) << '\t';
        short query_aa_index=q->get_first_aa_short_index(i-1);
        __m128i * sub_score_vec=(__m128i *)subMatrix->get_aa_substitution_vector(query_aa_index);
        __m128i score_matrix_vec01=*sub_score_vec;
        __m128i score_matrix_vec16=*(sub_score_vec+1);  
        size_t t_size=t->L;
        for (j=1; j<=t_size; ++j) // Loop through template positions j
            {
            
            // calculate amino acid profile-profile scores
            // __m128i score_vec=create_profile(&score_matrix_vec01,&score_matrix_vec16,&t->sequence[j-1]);
            
            
            __m128i template01=_mm_packs_epi16(t->sequence[j-1],t->sequence[j-1]);    
            // create slice mask
            // Example:
            //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
            //                      if lt 16 
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            __m128i lookup_mask01=_mm_cmplt_epi8(template01,sixteen_vec);
            __m128i lookup_mask16=_mm_cmpgt_epi8(template01,fiveteen_vec);
            // slice index
            // Example:
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
            //                          min	
            //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
            __m128i lookup_index01=_mm_min_epu8(lookup_mask01,template01);
            __m128i lookup_index16=_mm_min_epu8(lookup_mask16,template01);
            // 2xmal array lookup
            __m128i score01=_mm_shuffle_epi8(score_matrix_vec01,lookup_index01);
            __m128i score16 =_mm_shuffle_epi8(score_matrix_vec16,lookup_index16);
            
            __m128i score_vec_8bit=_mm_add_epi8(score01,score16);
            __m128i score_vec_16bit = _mm_srli_epi16(_mm_unpacklo_epi8(score_vec_8bit, score_vec_8bit), 8);
            
            
            
            
            score_vec_16bit = _mm_sub_epi16(score_vec_16bit,score_middle_vec);
            sDel_vec[i][j] = _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i][j-1],goe_vec),
                                           _mm_add_epi16(sDel_vec[i][j-1],ge_vec)
                                           );
            sIns_vec[i][j] = _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i-1][j],goe_vec),
                                           _mm_add_epi16(sIns_vec[i-1][j],ge_vec)
                                           );
            
            
            sMat_vec[i][j] = _mm_max_epi16(null_vec,
                                           _mm_max_epi16(
                                           _mm_add_epi16(sMat_vec[i-1][j-1],score_vec_16bit), 
                                           _mm_max_epi16(sDel_vec[i][j],sIns_vec[i][j])
                                                         )
                                           );
            col_max_score=_mm_max_epi16(col_max_score,sMat_vec[i][j]);
//            std::cout << (signed short) _mm_extract_epi16(sMat_vec[i][j],0) << '\t';

            
        } //end for j
//        std::cout << std::endl;	         

        // new score is higher
        // output 
        //  255 255 255 0   0   255 255 255
        // simulate Lower eq <=
        __m128i max_to_comp   = _mm_add_epi16(max_scores,minus_one_vec);
        __m128i lookup_mask_hi=_mm_cmplt_epi16(max_to_comp,col_max_score);
        // old score is higer
        // output 
        //  0   0   0   255 255   0   0   0
        __m128i lookup_mask_lo=_mm_cmpgt_epi16(max_scores,col_max_score);
        // -1  -1  -1   0   0  -1  -1  -1
        // -1  -1  -1  -1  -1  -1  -1  -1   
        //  1   1   1   0   0   1   1   1
        __m128i lookup_multp_hi=_mm_mullo_epi16(lookup_mask_hi,minus_one_vec);
        //  0   0   0  -1  -1   0   0   0
        // -1  -1  -1  -1  -1  -1  -1  -1   
        //  0   0   0   1   1   0   0   0
        __m128i lookup_multp_lo=_mm_mullo_epi16(lookup_mask_lo,minus_one_vec);
        //  1   1   1   0   0    1   1   1
        //                *
        //  s   s   s   0   0    s   s   s
        __m128i new_score_hi=_mm_mullo_epi16(lookup_multp_hi,col_max_score);
        
        //  0   0   0   1   1    0   0   0
        //                *
        //  0   0   0   s   s    0   0   0        
        __m128i new_score_lo=_mm_mullo_epi16(lookup_multp_lo,max_scores);
        //  0   0   0   s   s    0   0   0     
        //                +
        //  s   s   s   0   0    s   s   s
        //  s   s   s   s   s    s   s   s
        max_scores=_mm_add_epi16(new_score_hi,new_score_lo);
        __m128i curr_pos = _mm_set1_epi16(i);
        __m128i new_i_pos_hi=_mm_mullo_epi16(lookup_multp_hi,curr_pos);
        __m128i old_i_pos_lo=_mm_mullo_epi16(lookup_multp_lo,i_pos_max);
        i_pos_max= _mm_add_epi16(new_i_pos_hi,old_i_pos_lo);
        col_max_score=_mm_set1_epi16(NEGATIV_INF);
    } // end for i
    
    returnMatrix.scoreDel = this->scoreDel;
    returnMatrix.scoreMat = this->scoreMat;
    returnMatrix.scoreIns = this->scoreIns;
    short * max_scores_short=(short *) &max_scores;
    short * i_pos_max_short =(short *) &i_pos_max;
    for(size_t seq_index=0; seq_index<t->get_sequence_count(); seq_index++){
        size_t i=(size_t)i_pos_max_short[seq_index];
        returnMatrix.i[seq_index] = i;
        returnMatrix.j[seq_index] = 0;

        size_t j=t->get_sequence_size(seq_index);
        
        short max = (short) max_scores_short[seq_index];
        for(size_t col_index = j;col_index != 0;col_index--){
            short mat_col_score=scoreMat[i][col_index*(VEC_SIZE)+seq_index];
            if(max == mat_col_score){
                returnMatrix.j[seq_index] = col_index;
                break;
            }
        }
        
        returnMatrix.score[seq_index] =(float) (max /(float)subMatrix->scal_factor);
    }
    
    return returnMatrix;
}
////////////////////////////////////////////////////////////////////////




GotohEight::GotohMatrix GotohEight::global(SequenceEight* q, SequenceEight* t)
{
    const __m128i goe_vec= _mm_set1_epi16(goe);
    const __m128i ge_vec = _mm_set1_epi16(ge);
    const __m128i score_middle_vec=_mm_set1_epi16(subMatrix->middle_of_rang);
    __m128i ** sMat_vec = (__m128i **) (this->scoreMat);
    __m128i ** sDel_vec=(__m128i **)   (this->scoreDel);
    __m128i ** sIns_vec=(__m128i **)   (this->scoreIns);
    
    int i,j;      //query and template match state indices
    
    for (i=1; i<=q->L; ++i) // Loop through query positions i
    {
        
        // get score for query aa
        short query_aa_index=q->get_first_aa_short_index(i-1);
        __m128i * sub_score_vec=(__m128i *)subMatrix->get_aa_substitution_vector(query_aa_index);
        
        __m128i score_matrix_vec01=*sub_score_vec;
        __m128i score_matrix_vec16=*(sub_score_vec+1);  
        size_t t_size=t->L;
        for (j=1; j<=t_size; ++j) // Loop through template positions j
        {
            
            // calculate amino acid profile-profile scores
           // __m128i score_vec=create_profile(&score_matrix_vec01,&score_matrix_vec16,&t->sequence[j-1]);
            
 
            __m128i template01=_mm_packs_epi16(t->sequence[j-1],t->sequence[j-1]);    
            // create slice mask
            // Example:
            //	15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
            //                      if lt 16 
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            __m128i lookup_mask01=_mm_cmplt_epi8(template01,sixteen_vec);
            __m128i lookup_mask16=_mm_cmpgt_epi8(template01,fiveteen_vec);
            // slice index
            // Example:
            //  255	255	255	0	0	0	0	255	255	255	255	0	0	0	0	255
            //  15	12	11	16	20	19	18	11	15	12	11	16	20	19	18	11
            //                          min	
            //  15	12	11	0	0	0	0	11	15	12	11	0	0	0	0	155
            __m128i lookup_index01=_mm_min_epu8(lookup_mask01,template01);
            __m128i lookup_index16=_mm_min_epu8(lookup_mask16,template01);
            // 2xmal array lookup
            __m128i score01=_mm_shuffle_epi8(score_matrix_vec01,lookup_index01);
            __m128i score16 =_mm_shuffle_epi8(score_matrix_vec16,lookup_index16);
            
            __m128i score_vec_8bit=_mm_add_epi8(score01,score16);
            __m128i score_vec_16bit = _mm_srli_epi16(_mm_unpacklo_epi8(score_vec_8bit, score_vec_8bit), 8);
                
              
            
            
            score_vec_16bit = _mm_sub_epi16(score_vec_16bit,score_middle_vec);
            sDel_vec[i][j] = _mm_max_epi16(
                                     _mm_add_epi16(sMat_vec[i][j-1],goe_vec),
                                     _mm_add_epi16(sDel_vec[i][j-1],ge_vec)
                                     );
            sIns_vec[i][j] = _mm_max_epi16(
                                     _mm_add_epi16(sMat_vec[i-1][j],goe_vec),
                                     _mm_add_epi16(sIns_vec[i-1][j],ge_vec)
                                     );
            
            
            sMat_vec[i][j] = _mm_max_epi16(
                                     _mm_add_epi16(sMat_vec[i-1][j-1],score_vec_16bit), 
                                     _mm_max_epi16(sDel_vec[i][j],sIns_vec[i][j])
                                     );
            
//            std::cout << (signed short) _mm_extract_epi16(sMat_vec[i][j],0) << '\t';
 
            
        } //end for j
//        std::cout << std::endl;	         
        
    } // end for i
    
    returnMatrix.scoreDel = this->scoreDel;
    returnMatrix.scoreMat = this->scoreMat;
    returnMatrix.scoreIns = this->scoreIns;
    for(size_t seq_index=0; seq_index<t->get_sequence_count(); seq_index++){
        size_t i=q->get_sequence_size(seq_index);
        size_t j=t->get_sequence_size(seq_index);
        returnMatrix.i[seq_index] = i;
        returnMatrix.j[seq_index] = j;
        short score =scoreMat[i][j*(VEC_SIZE)+seq_index];
        returnMatrix.score[seq_index] =(float) (score /(float)subMatrix->scal_factor);
    }
    return returnMatrix;
}
////////////////////////////////////////////////////////////////////////

GotohEight::GotohMatrix GotohEight::calc_matrix(SequenceEight* q, SequenceEight* t){
    switch(ali_type) {
        case ALI_TYPE_LOCAL:
            return local(q, t);
        case ALI_TYPE_GLOBAL:
            return global(q, t);
        case ALI_TYPE_FREESHIFT:
        default:
            return freeshift(q, t);

    }
}



std::vector<char *> GotohEight::backtrace(GotohEight::GotohMatrix matrix,SequenceEight* q, SequenceEight* t){
    
    short ** score_del_matrix=matrix.scoreDel;
    short ** score_mat_matrix=matrix.scoreMat;
    short ** score_ins_matrix=matrix.scoreIns;
    std::vector<char *>  ret_vec;
    
    for(size_t seq_index = 0; seq_index < t->get_sequence_count(); seq_index++){
        size_t i = matrix.i[seq_index], j = matrix.j[seq_index];
        size_t current_index=0;

        
        size_t t_seq_size=t->get_sequence_size(seq_index);
        for(size_t seq_pos=t_seq_size; seq_pos!=j;seq_pos--){
            q_backtrace_seq_str[seq_index][current_index] = '-';
            t_backtrace_seq_str[seq_index][current_index] = t->get_char_by_pos(seq_pos-1, seq_index);
            ++current_index;
        }
    
        
        size_t q_seq_size=q->get_sequence_size(seq_index);
        for(size_t seq_pos=q_seq_size; seq_pos!=i;seq_pos--){
            t_backtrace_seq_str[seq_index][current_index] = '-';
            q_backtrace_seq_str[seq_index][current_index] = q->get_char_by_pos(seq_pos-1, seq_index);
            ++current_index;
        }
        
        
        
        while(i!=0|| j!=0){
            
            
            
            
            if(ali_type==ALI_TYPE_LOCAL){
                if(score_mat_matrix[i][(j)*(VEC_SIZE)+seq_index]==0){
                    break;
                }
            }
            
            if(i==0){
                q_backtrace_seq_str[seq_index][current_index] = '-';
                t_backtrace_seq_str[seq_index][current_index] = t->get_char_by_pos(j-1, seq_index);
                current_index++;
                j--;
            }else 
            if(j==0){
                q_backtrace_seq_str[seq_index][current_index] = q->get_char_by_pos(i-1, seq_index);
                t_backtrace_seq_str[seq_index][current_index] = '-';
                current_index++;
                i--;
            } else { 
                // Ai,j = Ai-1,j-1 + S(si,tj) 
            const short t_aa_index=t->get_aa_short_index(j-1, seq_index);
            const short q_aa_index=q->get_first_aa_short_index(i-1);
            const short ss_score=subMatrix->get_substituion_score(q_aa_index,t_aa_index)-subMatrix->middle_of_rang;
            if(score_mat_matrix[i-1][(j-1)*VEC_SIZE+seq_index]+ss_score==score_mat_matrix[i][(j)*(VEC_SIZE)+seq_index]){
                    q_backtrace_seq_str[seq_index][current_index] = q->get_char_by_pos(i-1, seq_index);
                    t_backtrace_seq_str[seq_index][current_index] = t->get_char_by_pos(j-1, seq_index);
                    current_index++;
                    i--;
                    j--;
                    
            }else   
                // Ai,j = Ii,j
                if(score_mat_matrix[i][(j)*VEC_SIZE+seq_index] == score_ins_matrix[i][(j)*VEC_SIZE+seq_index]){
                    int counter = 0;
                    do{
                        counter++;
                        q_backtrace_seq_str[seq_index][current_index] = q->get_char_by_pos(i-counter, seq_index);
                        t_backtrace_seq_str[seq_index][current_index] = '-';
                        current_index++;
                    }while(score_mat_matrix[i][(j)*VEC_SIZE+seq_index]
                           !=
                           score_mat_matrix[i-counter][(j)*VEC_SIZE+seq_index]+goe+(counter-1)*ge );   
                    i=i-counter;
                }else  
                    // Ai,j = Di,j 
                    if(score_mat_matrix[i][(j)*VEC_SIZE+seq_index] == score_del_matrix[i][(j)*VEC_SIZE+seq_index]){
                        int counter = 0;
                        do{
                            counter++;
                            q_backtrace_seq_str[seq_index][current_index] = '-';
                            t_backtrace_seq_str[seq_index][current_index] = t->get_char_by_pos(j-counter, seq_index);
                            current_index++;
                        }while(score_mat_matrix[i][(j)*VEC_SIZE+seq_index]
                               !=
                               score_mat_matrix[i][(j-counter)*VEC_SIZE+seq_index]+goe+(counter-1)*ge );
                        j=j-counter;
                    }

            }         
        }
        for(size_t seq_pos=j; seq_pos!=0;seq_pos--){
            q_backtrace_seq_str[seq_index][current_index] = '-';
            t_backtrace_seq_str[seq_index][current_index] = t->get_char_by_pos(seq_pos-1, seq_index);
            ++current_index;
        }
        
        for(size_t seq_pos=i; seq_pos!=0;seq_pos--){
            t_backtrace_seq_str[seq_index][current_index] = '-';
            q_backtrace_seq_str[seq_index][current_index] = q->get_char_by_pos(seq_pos-1, seq_index);
            ++current_index;
        }
        
        q_backtrace_seq_str[seq_index][current_index]='\0';
        t_backtrace_seq_str[seq_index][current_index]='\0';
        strrev(q_backtrace_seq_str[seq_index]);
        strrev(t_backtrace_seq_str[seq_index]);
        ret_vec.push_back(q_backtrace_seq_str[seq_index]);
        ret_vec.push_back(t_backtrace_seq_str[seq_index]);
                
    }
    
    return (ret_vec);
}

float GotohEight::calc_check_score(char *q, char * t) {

    size_t start = 0;
	size_t end = strlen(q);
	if(this->ali_type ==ALI_TYPE_LOCAL || ali_type == ALI_TYPE_FREESHIFT) {
		bool found_first = false;
		for(size_t i = 0; i < strlen(q); i++) {
            if(q[i] != '-' && t[i] != '-') {
                if(!found_first) {
					start = i;
					found_first = true;
				}
				end = i + 1;
			}
		}
	}
    
    int score = 0;
    bool is_gap_open = false;
    for(size_t i = start; i < end; i++) {
        if(q[i] == '-' || t[i] == '-') {
            if(is_gap_open) {
                score += this->ge;
            } else {
                score += this->goe;
                is_gap_open = true;
            }
            continue;
        }
        score += subMatrix->get_substituion_score(q[i] , t[i] )-subMatrix->middle_of_rang;
        is_gap_open = false;
    }
    
    return ((float)score) / ((float) subMatrix->scal_factor);
}

