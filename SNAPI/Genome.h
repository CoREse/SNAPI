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
		unsigned long long start_location;
		NASeq sequence;
	}
	Genome();
	~Genome();
	std::vector<Chromesome> chrs;//the start_locations of those chromesomes are increasing, that's to say, the earlier you are read, the smaller your start_location is
	unsigned getChrNumber(unsigned long long location) const
	{
		for (unsigned i=0;i<MAX_CHROMESOME;++i)
		{
			if (location<end_locations) return i;
		}
	}
private:
	unsigned long long start_locations[MAX_CHROMESOMES];
};

