#pragma once
#include "NASeq.h"
#include <string>
#include <vector>

#define MAX_CHROMESOMES 100

class Genome
{
public:
	struct Chromesome
	{
		std::string name;
		unsigned start_location;
		NASeq sequence;
	};
	Genome(const char * FolderPath);//take the folder of fa files as the parameter
	~Genome();
	std::vector<Chromesome> chrs;//the start_locations of those chromesomes are increasing, that's to say, the earlier you are read, the smaller your start_location is
	unsigned getChrNumber(unsigned location) const
	{
		for (unsigned i = 0; i < MAX_CHROMESOMES; ++i)
		{
			if (location < end_locations[i]) return i;
		}
	}
	bool saveToFile(FILE*);
	bool readFromFile(FILE*);
private:
	unsigned end_locations[MAX_CHROMESOMES];
	friend class Indexer;
};

