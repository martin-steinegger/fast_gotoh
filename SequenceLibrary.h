#pragma once
#include <string>
#include <map>

class SequenceLibrary
{
public:
	SequenceLibrary(std::string file);
	std::pair<std::string,size_t> get_sequence(std::string id);
    size_t get_max_seq_size();
private:
	std::map<std::string, std::pair< std::string,size_t> > sequences;
    size_t max_seq_size;
};

