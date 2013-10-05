    //
//  main.cpp
//  run
//
//  Created by Martin Steinegger on 23.10.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include <iostream>
#include "ParseMatrix.h"
#include "SubstitutionMatrix.h"
#include "SequenceEight.h"
#include "GotohEight.h"
#include "PairsLibrary.h"
#include "SequenceLibrary.h"
#include <algorithm>
#include "util.h"
#include <stdio.h>
#include <string.h>

bool sort_seq_vector (std::pair<std::string, int> i,std::pair<std::string, int> j);
bool sort_seq_vector (std::pair<std::string, int> i,std::pair<std::string, int> j) { return (i.second > j.second); }

inline void print_score(GotohEight::GotohMatrix matrix,const std::string query_key,const char** template_keys,int elements);
inline void print_score(GotohEight::GotohMatrix matrix,const std::string query_key,const char** template_keys,int elements){

    for(size_t i = 0; i < elements; i++) {
        printf("%s %s %.3f\n",query_key.c_str(),template_keys[i],matrix.score[i]);
    }

}


inline void print_ali_result(char * q_str,char* t_str,float score,
                      const std::string query_key,const char* template_key);
inline void print_ali_result(char * q_str,char* t_str,float score,
                      const std::string query_key,const char* template_key){
    printf(">%s %s %.3f\n",query_key.c_str(),template_key,score);
    std::cout << query_key <<": "<< q_str <<'\n';
    std::cout << template_key <<": "<< t_str <<'\n';
}

void help(char all);
void help(char all=0)
{
    printf("\n");
    printf("Parameter:\n");
    printf("  -seqlib         path to sequence libary\n");
    printf("  -pairs          path to pairs file\n");
    printf("  -m              matrixname (Default dayhoff)\n");
    printf("  -go             gapopen (Default -12) as positiv number\n");
    printf("  -ge             gapextend (Default -1) as positiv number\n");
    printf("  -mode           local|global|freeshift (Default freeshift)\n");
    printf("  -printali       print the Alignment\n");
    printf("  -check          checks Scores based on the alignments(will print ali to)\n");
}
struct ProgramArgument {
    std::string matrixname;
    std::string seqlib;
    std::string pairs;
    int go;
    int ge;
    int mode;
    bool print_ali;
    bool printmatrices;
    bool check;
};

/////////////////////////////////////////////////////////////////////////////////////
//// Processing input options from command line 
/////////////////////////////////////////////////////////////////////////////////////
ProgramArgument ProcessArguments(int argc,const char** argv);
ProgramArgument ProcessArguments(int argc,const char** argv)
{
    ProgramArgument program_args;
    
    program_args.print_ali = false;
    program_args.printmatrices = false;
    program_args.mode = ALI_TYPE_FREESHIFT;
    program_args.go = -12;
    program_args.ge = -1;
    program_args.matrixname = "matrices/dayhoff.mat";
    program_args.check = false;
    //Processing command line input
    if (argc == 1){
        help();
        exit(4);
    }
    for (int i=1; i<argc; i++)
    {
        if (!strcmp(argv[i],"-pairs"))
        {
            if (++i>=argc || argv[i][0]=='-')
            {help() ; std::cerr<< std::endl <<"Error in gotoh: no pairs file  -pairs\n"; exit(4);}
            else 
                program_args.pairs = argv[i];
        }
        else if (!strcmp(argv[i],"-seqlib"))
        {
            if (++i>=argc || argv[i][0]=='-')
            {help() ; std::cerr<< std::endl <<"Error in gotoh: no sequence libary -seqlib\n"; exit(4);}
            else 
                program_args.seqlib = argv[i];
        }
        else if (!strcmp(argv[i],"-go"))
        {
            if (++i>=argc || argv[i][0]=='-')
                {help() ; std::cerr<< std::endl <<"Error in gotoh: no gap open costs -go\n"; exit(4);}
            else 
                program_args.go = -atoi(argv[i]);
        }
        else if (!strcmp(argv[i],"-ge"))
        {
            if (++i>=argc || argv[i][0]=='-')
                {help() ; std::cerr<< std::endl <<"Error in gotoh: no gap extend costs -ge\n"; exit(4);}
            else
                program_args.ge = -atoi(argv[i]);
        }
        else if (!strcmp(argv[i],"-printali"))
        {
            program_args.print_ali = true;
        }
        else if (!strcmp(argv[i],"-m"))
        {
            if (++i>=argc || argv[i][0]=='-')
                {help() ; std::cerr<< std::endl <<"Error in gotoh: no matrix name -m\n"; exit(4);}
            else
                program_args.matrixname = std::string(argv[i]);
        }
        else if (!strcmp(argv[i],"-mode"))
        {
            if (++i>=argc || argv[i][0]=='-')
            {help() ; std::cerr<< std::endl <<"Error in gotoh: no mode name -mode\n"; exit(4);}
            else{
                std::string mode(argv[i]);
                if(mode.compare("freeshift")==0){
                    program_args.mode = ALI_TYPE_FREESHIFT;
                }else if(mode.compare("global")==0){
                    program_args.mode = ALI_TYPE_GLOBAL;
                }else if(mode.compare("local")==0){
                    program_args.mode = ALI_TYPE_LOCAL;
                }
            }
        }

        
        else if (!strcmp(argv[i],"-check"))
        {
            
            program_args.check = true;
            program_args.print_ali = true;
        }
        else if (!strcmp(argv[i],"-printmatrices"))
        {
            
            program_args.printmatrices = true;
        }
        
    }
    return program_args;
}


int main (int argc, const char ** argv)
{

    
    ProgramArgument program_args=ProcessArguments(argc,argv);
    
    
    bool print_ali=program_args.print_ali;

    ParseMatrix submatParser(program_args.matrixname);
#ifdef DEBUG
    std::cout << "Print Matrix\n";
    submatParser.print();
#endif
    SubstitutionMatrix subMatrix(submatParser.scores, submatParser.row_index);
#ifdef DEBUG
    std::cout << "Print SubMatrix!\n";
    subMatrix.print();
#endif
    SequenceLibrary sequences(program_args.seqlib);
    const size_t max_seq_size = sequences.get_max_seq_size();
    PairsLibrary pairs(program_args.pairs,&sequences);
    GotohEight gotohAlgo(max_seq_size,&subMatrix, program_args.go, program_args.ge, program_args.mode);
    SequenceEight sequenceQuery(max_seq_size, subMatrix.aa2short,subMatrix.short2aa);
    SequenceEight sequenceTemplate(max_seq_size, subMatrix.aa2short,subMatrix.short2aa);

    
        
    char buffer[2097152];
    std::cout.rdbuf()->pubsetbuf(buffer, 2097152);
    
    const char ** template_sequences=(const char **)memalign(16,(VEC_SIZE)*sizeof(const char *)); 
    const char ** template_keys=(const char **)memalign(16,(VEC_SIZE)*sizeof(const char *)); 
    
    std::map< std::string,std::vector<std::pair<std::string,size_t> > >::iterator it = pairs.pairs.begin();    
    for( ; it != pairs.pairs.end(); ++it ){
        std::pair<std::string, size_t> query = sequences.get_sequence(it->first);
        std::string query_key = it->first.c_str();
        sequenceQuery.MapOneSEQToSEQ8(query.first.c_str());
        std::vector<std::pair<std::string, size_t> > t_seq = it->second;
        sort (t_seq.begin(), t_seq.end(), sort_seq_vector); 

        int elem_count=0;
        size_t t_seq_size = t_seq.size();
        for(std::vector<int>::size_type i = 0; i < t_seq_size ; i++) {
            std::string template_key=t_seq[i].first;
            template_keys[i%VEC_SIZE] = template_key.c_str();
            template_sequences[i%VEC_SIZE] = sequences.get_sequence(template_key).first.c_str();
            elem_count++;
            if((i+1)%VEC_SIZE==0){
                elem_count=0;
                sequenceTemplate.MapSEQArrayToSEQ8(template_sequences, VEC_SIZE);
                GotohEight::GotohMatrix matrix= gotohAlgo.calc_matrix(&sequenceQuery, &sequenceTemplate);
                if(print_ali == true){
                    std::vector<char *> alignments=gotohAlgo.backtrace(matrix,&sequenceQuery,&sequenceTemplate);
                        for(size_t i = 0; i < 8; i++) {
                            if(program_args.check==true){
                                float check_score=gotohAlgo.calc_check_score(alignments[(i*2)],alignments[(i*2)+1]);
                                printf("check score: %.3f\n",check_score);
                            }
                            print_ali_result(alignments[(i*2)],alignments[(i*2)+1],matrix.score[i],query_key,template_keys[i]);
                        }
                }else{    
                    print_score(matrix,query_key,template_keys,VEC_SIZE);
                }
            }
            
        }
        if(elem_count!=0){
            sequenceTemplate.MapSEQArrayToSEQ8(template_sequences, elem_count);
            GotohEight::GotohMatrix matrix= gotohAlgo.calc_matrix(&sequenceQuery, &sequenceTemplate);
            if(print_ali == true){
                std::vector<char *> alignments=gotohAlgo.backtrace(matrix,&sequenceQuery,&sequenceTemplate);
                for(size_t i = 0; i < elem_count; i++) {
                    if(program_args.check==true){
                        float check_score=gotohAlgo.calc_check_score(alignments[(i*2)],alignments[(i*2)+1]);
                        printf("check score: %.3f\n",check_score);
                    }
                print_ali_result(alignments[(i*2)],alignments[(i*2)+1],matrix.score[i],query_key,template_keys[i]);
                }
            }else{    
                print_score(matrix,query_key,template_keys,elem_count);
            }
        }
        
    }
    std::cout.flush();
    free(template_keys);
    free(template_sequences);
    
    
    
    
    return 0;
}

