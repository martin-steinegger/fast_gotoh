#pragma once
#include <string>
#include <vector>
#include <map>
#include "SequenceLibrary.h"

class PairsLibrary
{
public:
	PairsLibrary(std::string file,SequenceLibrary * seq_lib);
    
    
    std::map< std::string,std::vector<std::pair<std::string,size_t> > > pairs;

};

