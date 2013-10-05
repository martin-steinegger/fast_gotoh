#include "PairsLibrary.h"
#include <fstream>
#include "util.h"


PairsLibrary::PairsLibrary(std::string file,SequenceLibrary * seq_lib)
{
	std::ifstream in(file.c_str());
	std::string line;
	while(std::getline(in, line)) {
        std::vector<std::string> myvect= split(line,' ');
        std::string key = myvect.at(0);
        size_t seq_size=seq_lib->get_sequence(key).second;
        std::pair<std::string,size_t> to_add(myvect.at(1),seq_size);
        
        
        pairs[key].push_back(to_add);
	}
}
