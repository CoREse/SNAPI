#include "Indexer.h"
#include "seed.h"

Genome * Indexer::Reference = nullptr;
HashTable* Indexer::Index = nullptr;
unsigned Indexer::MaxH = 0;

Indexer::Indexer()
{
}


Indexer::~Indexer()
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
	fprintf(file, "%u\n", Reference->chrs.size());
	//first the reference's chrs
	for (unsigned i = 0; i < Reference->chrs.size(); ++i)
	{
		fprintf(file, "%s %u\n", Reference->chrs[i].name, Reference->chrs[i].start_location);
		Reference->chrs[i].sequence.saveToFile(file);
		fprintf(file, "\n");
		fprintf(file,"%u\n",Reference->end_locations[i]);
	}
	fclose(file);
	
	//create the table filenames
	const char *tableFileExtension = ".table";
	size_t tableFileNameLength = (size_t)(strlen(fname) + strlen(tableFileExtension) + 4);
	char *tableFileName = new char[tableFileNameLength];
	int No = 0;
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
			throw - 110;//–¥»ÎIndex ß∞‹£°
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
	fprintf(file, "%u\n", Index->nOverflowTable);
	//first the reference's chrs
	for (unsigned i = 0; i < Index->nOverflowTable; ++i)
	{
		fprintf(file, "%ull\n", Index->OverflowTable[i].key);
		for (unsigned j = 0; j < Index->OverflowTable[i].locations.size(); ++j)
		{
			fprintf(file, "%u", Index->OverflowTable[i].locations[j]);
		}
	}
	fclose(file);
}

bool Indexer::loadFromFile(const char * fname)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file = fopen(fname, "r");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	fscanf(file, "%s%u%u%u%u", name, &size, &lengthOfGenome, &usedIndex, &seed::length);
	fclose(file);
	//allocate space
	free(Index);
	Index = (entry*)malloc(size*sizeof(entry));
	if (Index == nullptr)
	{
		fclose(file);
		throw - 106;//malloc∑÷≈‰ø’º‰ ß∞‹
		return false;
	}
	//create the table filenames
	const char *tableFileExtension = ".table";
	size_t tableFileNameLength = (size_t)(strlen(fname) + strlen(tableFileExtension) + 4);
	char *tableFileName = new char[tableFileNameLength];
	int No = 0;
	do
	{
		_snprintf(tableFileName, tableFileNameLength, "%s%s%03d", fname, tableFileExtension, No);
		//read Index from the table files
		file = fopen(tableFileName, "rb");
		if (file == nullptr)
		{
			throw - 102;
			return false;
		}
		if (fread(Index + MAX_SAVE_SIZE*No, sizeof(entry), size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE, file) != (size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE))
		{
			fclose(file);
			throw - 111;//∂¡»°Index ß∞‹£°
			return false;
		}
		fclose(file);
		++No;
		if (size <= MAX_SAVE_SIZE*No) break;
	} while (No < 1000);
	//create the genome filename
	const char *genomeFileExtension = ".genome";
	size_t genomeFileNameLength = (size_t)(strlen(fname) + strlen(genomeFileExtension) + 1);
	char *genomeFileName = new char[genomeFileNameLength];
	_snprintf(genomeFileName, genomeFileNameLength, "%s%s", fname, genomeFileExtension);
	//allocate space
	free(theGenome);
	Reference = (char*)malloc(lengthOfGenome*sizeof(char));
	if (theGenome == nullptr)
	{
		fclose(file);
		throw - 106;//malloc∑÷≈‰ø’º‰ ß∞‹
		return false;
	}
	//read theGenome from the genome file
	file = fopen(genomeFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	Reference = new Genome();
	unsigned Gsize;
	fscanf(file, "%u", &Gsize);
	extern char * buffer;
	for (unsigned i = 0; i < Gsize; ++i)
	{
		Genome::Chromesome tmpChr;
		fscanf(file, "%s %u\n", buffer, tmpChr.start_location);
		tmpChr.name = buffer;
		tmpChr.sequence.loadFromFile(file);
		Reference->chrs.push_back(tmpChr);
		fprintf(file, "\n");
	}
	fclose(file);
	if (fread(theGenome, sizeof(char), lengthOfGenome, file) != lengthOfGenome)
	{
		fclose(file);
		throw - 111;//∂¡»°Index ß∞‹£°
		return false;
	}
	fclose(file);
	return true;
}