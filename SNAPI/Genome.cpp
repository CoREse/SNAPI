#include "Genome.h"
#include <io.h>
#include <stdio.h>
#include <cassert>

using namespace std;

char Genome::NameBuffer[1024];

Genome::Genome(const char * FolderPath)
{
	_finddata_t file;
	intptr_t lf;
	string FindString = FolderPath;
	FindString.append("*.fa");
	if ((lf = _findfirst(FindString.c_str(), &file)) != -1l)
	{
		FILE* infile;
		Chromesome tmpChr;
		do
		{
			assert(chrs.size() < MAX_CHROMESOMES);
			//open the fa file
			string FaName = FolderPath;
			FaName.append(file.name);
			infile = fopen(FaName.c_str(),"r");
			if (infile == nullptr)
			{
				throw - 102;
				return;
			}
			//read the fa file
			fscanf(infile, "%s", NameBuffer);
			tmpChr.name = NameBuffer;
			rewind(infile);
			NASeq tmpNASeq(infile);
			fclose(infile);
			//push back a new chromesome
			tmpChr.sequence = tmpNASeq;
			if (chrs.size() == 0)
			{
				tmpChr.start_location = 0;
			}
			else
			{
				tmpChr.start_location = end_locations[chrs.size() - 1];
			}
			end_locations[chrs.size()] = tmpChr.start_location + tmpChr.sequence.length;
			chrs.push_back(tmpChr);
		} while (_findnext(lf,&file)==0);
	}
	_findclose(lf);
}

Genome::Genome()
{
}

Genome::~Genome()
{
}

bool Genome::saveToFile(FILE* file)
{
	fprintf(file,"%u\n",chrs.size());
	for (unsigned i=0;i<chrs.size();++i)
	{
		fprintf(file,"%s ",chrs[i].name.c_str());
		if (!chrs[i].sequence.saveToFile(file)) return false;
		fprintf(file,"%u %u\n",chrs[i].start_location, end_locations[i]);
	}
	return true;
}

bool Genome::loadFromFile(FILE* file)
{
	unsigned chrsize;
	fscanf(file,"%u\n",&chrsize);
	for (unsigned i=0;i<chrsize;++i)
	{
		Chromesome tmpChr;
		fscanf(file,"%s ",NameBuffer);
		tmpChr.name=NameBuffer;
		if (!tmpChr.sequence.loadFromFile(file)) return false;
		fscanf(file,"%u %u\n",&(tmpChr.start_location),end_locations+i);
		chrs.push_back(tmpChr);
	}
	return true;
}
