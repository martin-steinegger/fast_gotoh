#include "SequenceLibrary.h"
#include <fstream>
#include "util.h"


SequenceLibrary::SequenceLibrary(std::string file)
{
	std::ifstream in(file.c_str());
	std::string line;
    max_seq_size = 0;
	while(std::getline(in, line)) {
        std::vector<std::string> myvect= split(line,':');
        size_t seq_size=myvect.at(1).length();
        max_seq_size=std::max(seq_size, max_seq_size);
        std::string key = myvect.at(0);
        std::pair<std::string,size_t> to_add( myvect.at(1),seq_size);
		sequences.insert(std::pair<std::string, std::pair<std::string,int> >
                             (key,to_add));
	}
}

std::pair<std::string,size_t> SequenceLibrary::get_sequence(std::string id) {
	return sequences.find(id)->second;
}

size_t SequenceLibrary::get_max_seq_size(){
    return max_seq_size;
}

