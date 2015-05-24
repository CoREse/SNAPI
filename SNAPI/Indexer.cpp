#include "Indexer.h"
#include "seed.h"
#include <cassert>
#include "BigAlloc.h"

Genome * Indexer::Reference = nullptr;
HashTable* Indexer::Index = nullptr;
unsigned Indexer::MaxH = 0;
unsigned Indexer::ReadLength = 100;

Indexer::Indexer()
{
}


Indexer::~Indexer()
{
}

void Indexer::cleanIndex()
{
	if (Reference != nullptr) delete Reference;
	if (Index != nullptr) delete Index;
}

bool Indexer::saveToFile(const char * fname)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file;
	//create the genome filename
	const char *genomeFileExtension = ".genome";
	size_t genomeFileNameLength = (size_t)(strlen(fname) + strlen(genomeFileExtension) + 1);
	char *genomeFileName = new char[genomeFileNameLength];
	_snprintf(genomeFileName, genomeFileNameLength, "%s%s", fname, genomeFileExtension);
	//save the Index to the genome file
	file = fopen(genomeFileName, "wb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	if(!Reference->saveToFile(file)) return false;
	fclose(file);
	
	//create the table filenames
	const char *tableFileExtension = ".table";
	size_t tableFileNameLength = (size_t)(strlen(fname) + strlen(tableFileExtension) + 4);
	char *tableFileName = new char[tableFileNameLength];
	int No = 0;
	_snprintf(tableFileName, tableFileNameLength, "%s%s", fname, tableFileExtension);
	//save the Index to the table files
	file = fopen(tableFileName, "wb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	fprintf(file, "%f %f\n",HashTable::MainPara, HashTable::OverflowPara);
	fprintf(file, "%u %u %u %u\n", Index->nMainTable, Index->nUsedMainTable, Index->nOverflowTable, Index->nUsedOverflowTable);
	do
	{
		_snprintf(tableFileName, tableFileNameLength, "%s%s%03d", fname, tableFileExtension, No);
		//save the Index to the table files
		file = fopen(tableFileName, "wb");
		if (file == nullptr)
		{
			throw - 102;
			return false;
		}
		if (fwrite(Index->MainTable + MAX_SAVE_SIZE*No, sizeof(HashTable::Entry), Index->nMainTable <= MAX_SAVE_SIZE*(No + 1) ? Index->nMainTable - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE, file) != (Index->nMainTable <= MAX_SAVE_SIZE*(No + 1) ? Index->nMainTable - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE))
		{
			fclose(file);
			throw - 110;//Ð´ÈëIndexÊ§°Ü£¡
			return false;
		}
		fclose(file);
		++No;
		if (Index->nMainTable <= MAX_SAVE_SIZE*No) break;
	} while (No < 1000);
	//create the overflow table filenames
	const char *oftFileExtension = ".overflowtable";
	size_t oftFileNameLength = (size_t)(strlen(fname) + strlen(oftFileExtension) + 1);
	char *oftFileName = new char[oftFileNameLength];
	_snprintf(oftFileName, oftFileNameLength, "%s%s", fname, oftFileExtension);
	//save the Index to the overflow table
	file = fopen(oftFileName, "wb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	for (unsigned i = 0; i<Index->nOverflowTable; ++i)
	{
		fprintf(file, "%llu %u ", Index->OverflowTable[i].key, Index->OverflowTable[i].locations.size());
		for (unsigned j = 0; j<Index->OverflowTable[i].locations.size(); ++j)
		{
			fprintf(file, "%u ", Index->OverflowTable[i].locations[j]);
		}
	}
	fclose(file);
	return true;
}

bool Indexer::loadFromFile(const char * fname)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file;
	const char *tableFileExtension = ".table";
	size_t tableFileNameLength = (size_t)(strlen(fname) + strlen(tableFileExtension) + 4);
	char *tableFileName = new char[tableFileNameLength];
	int No = 0;
	_snprintf(tableFileName, tableFileNameLength, "%s%s", fname, tableFileExtension);
	//save the Index to the table files
	file = fopen(tableFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	Index = new HashTable();
	fscanf(file, "%f %f\n", &(HashTable::MainPara), &(HashTable::OverflowPara));
	fscanf(file, "%u %u %u %u\n", &(Index->nMainTable), &(Index->nUsedMainTable), &(Index->nOverflowTable), &(Index->nUsedOverflowTable));
	fclose(file);
	Index->MainTable = (HashTable::Entry*)malloc(Index->nMainTable*sizeof(HashTable::Entry));
	do
	{
		_snprintf(tableFileName, tableFileNameLength, "%s%s%03d", fname, tableFileExtension, No);
		//save the Index to the table files
		file = fopen(tableFileName, "rb");
		if (file == nullptr)
		{
			throw - 102;
			return false;
		}
		if (fread(Index->MainTable + MAX_SAVE_SIZE*No, sizeof(HashTable::Entry), Index->nMainTable <= MAX_SAVE_SIZE*(No + 1) ? Index->nMainTable - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE, file) != (Index->nMainTable <= MAX_SAVE_SIZE*(No + 1) ? Index->nMainTable - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE))
		{
			fclose(file);
			throw - 110;//Ð´ÈëIndexÊ§°Ü£¡
			return false;
		}
		fclose(file);
		++No;
		if (Index->nMainTable <= MAX_SAVE_SIZE*No) break;
	} while (No < 1000);
	//read the overflow Table
	const char *oftFileExtension = ".overflowtable";
	size_t oftFileNameLength = (size_t)(strlen(fname) + strlen(oftFileExtension) + 1);
	char *oftFileName = new char[oftFileNameLength];
	_snprintf(oftFileName, oftFileNameLength, "%s%s", fname, oftFileExtension);
	file = fopen(oftFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	unsigned nloc, tmploc;
	Index->OverflowTable = (HashTable::OverflowEntry*)BigAlloc(Index->nOverflowTable*sizeof(HashTable::OverflowEntry));
	for (unsigned i = 0; i<Index->nOverflowTable; ++i)
	{
		fscanf(file, "%llu %u ", &(Index->OverflowTable[i].key), &nloc);
		for (unsigned j = 0; j<nloc; ++j)
		{
			fscanf(file, "%u ", &tmploc);
			Index->OverflowTable[i].locations.push_back(tmploc);
		}
	}
	fclose(file);

	//create the genome filename
	const char *genomeFileExtension = ".genome";
	size_t genomeFileNameLength = (size_t)(strlen(fname) + strlen(genomeFileExtension) + 1);
	char *genomeFileName = new char[genomeFileNameLength];
	_snprintf(genomeFileName, genomeFileNameLength, "%s%s", fname, genomeFileExtension);
	assert(Reference==nullptr);
	file = fopen(genomeFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	Reference = new Genome();
	Reference->loadFromFile(file);
	fclose(file);
	return true;
}

Genome* Indexer::readReference(const char * FolderPath)
{
	assert(Reference == nullptr);
	Reference = new Genome(FolderPath);
	return Reference;
}

HashTable* Indexer::buildIndex()
{
	if (Index != nullptr) delete Index;
	Index = new HashTable(Reference->getTotalLength());
	for (unsigned i = 0; i != Reference->chrs.size(); ++i)
	{
		for (unsigned j = 0; j != Reference->chrs[i].sequence.length; ++j)
		{
			NASeq tmpSeq = Reference->chrs[i].sequence.getSubSequence(j, seed::seedLen);
			if (seed::isASeed(tmpSeq))
			{
				seed TheSeed(tmpSeq);
				Index->insert(TheSeed.getBases(), Reference->chrs[i].start_location + j);
			}
		}
	}
	return Index;
}

HashTable* Indexer::cutIndex(unsigned hMax)
{
	HashTable::Entry * tmpEntry;
	unsigned tmpLocation[2],tmpStartLoc[2],tmpChrNum[2];
	int tmpEdd;
	for (unsigned i = 0; i != Index->nOverflowTable; ++i)
	{
		if (Index->OverflowTable[i].locations.size() > hMax - 2)
		{
			tmpEntry = Index->lookupMainTable(Index->OverflowTable[i].key);
			tmpLocation[0] = tmpEntry->location1;
			for (auto j = Index->OverflowTable->locations.begin(); j != Index->OverflowTable->locations.end(); ++j)
			{
				tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
				tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
				tmpLocation[1] = Index->OverflowTable->locations[*j];
				tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
				tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
				tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.getSubSequence(tmpLocation[0] - tmpStartLoc[0] - ReadLength, 2 * ReadLength),
					Reference->chrs[tmpChrNum[1]].sequence.getSubSequence(tmpLocation[1] - tmpStartLoc[1] - ReadLength, 2 * ReadLength),
					MAX_EDD
					);
				if (tmpEdd > 0 && tmpEdd < 10)//consider as similiar
				{
					Index->OverflowTable->locations.erase(j);
				}
			}
			tmpLocation[0] = tmpEntry->location1;
			for (auto j = Index->OverflowTable->locations.begin(); j != Index->OverflowTable->locations.end(); ++j)
			{
				tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
				tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
				tmpLocation[1] = Index->OverflowTable->locations[*j];
				tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
				tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
				tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.getSubSequence(tmpLocation[0] - tmpStartLoc[0] - ReadLength, 2 * ReadLength),
					Reference->chrs[tmpChrNum[1]].sequence.getSubSequence(tmpLocation[1] - tmpStartLoc[1] - ReadLength, 2 * ReadLength),
					MAX_EDD
					);
				if (tmpEdd > 0 && tmpEdd < 10)//consider as similiar
				{
					Index->OverflowTable->locations.erase(j);
				}
			}
			for (auto q = Index->OverflowTable->locations.begin(); q != Index->OverflowTable->locations.end(); ++q)
			{
				tmpLocation[0] = *q;
				for (auto j = Index->OverflowTable->locations.begin(); j != Index->OverflowTable->locations.end(); ++j)
				{
					tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
					tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
					tmpLocation[1] = Index->OverflowTable->locations[*j];
					tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
					tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
					tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.getSubSequence(tmpLocation[0] - tmpStartLoc[0] - ReadLength, 2 * ReadLength),
						Reference->chrs[tmpChrNum[1]].sequence.getSubSequence(tmpLocation[1] - tmpStartLoc[1] - ReadLength, 2 * ReadLength),
						MAX_EDD
						);
					if (tmpEdd > 0 && tmpEdd < 10)//consider as similiar
					{
						Index->OverflowTable->locations.erase(j);
					}
				}
			}
		}
	}
	return Index;
}