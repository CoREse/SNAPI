#include "Genome.h"
#include <io.h>
#include <stdio.h>
#include <cassert>

using namespace std;

Genome::Genome(const char * FolderPath)
{
	_finddata_t file;
	long lf;
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
			extern char *buffer;//defined in NASeq.cpp
			string FaName = FolderPath;
			FaName.append(file.name);
			infile = fopen(FaName.c_str(),"r");
			if (infile == nullptr)
			{
				throw - 102;
				return;
			}
			//read the fa file
			fscanf(infile, "%s", buffer);
			tmpChr.name = buffer;
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


Genome::~Genome()
{
}
