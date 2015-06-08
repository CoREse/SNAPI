#include "Indexer.h"
#include "seed.h"
#include <cassert>
#include "BigAlloc.h"
#include <string.h>
#include <direct.h>
#include <iostream>
#include <omp.h>
#include "Aligner.h"

Genome * Indexer::Reference = nullptr;
HashTable* Indexer::Index = nullptr;
unsigned Indexer::hMax = 100;
unsigned Indexer::ReadLength = 100;
unsigned Indexer::SimilarEdd = 10;

#define TO_BE_CLEAR 2//if we found TO_BE_CLEAR number of similar locations, we remove all of them instead leave one reserved
#define PERSAVE 500000

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

bool Indexer::saveToFile(const char * Fname)
{
	if (Fname == nullptr)
	{
		return false;
	}
	FILE * file;
	char* slashI;
	char fname[MAX_NAME_LENGTH];
	strcpy(fname, Fname);
#ifdef WINDOWS
	slashI = strrchr(fname, '\\');
	if (slashI != nullptr)
	{
		*slashI = '\0';
		_mkdir(fname);
		*slashI = '\\';
	}
#endif
#ifdef LINUX
	slashI = strrchr(fname, '/');
	if (slashI != nullptr)
	{
		*slashI = '\0';
		_mkdir(fname);
		*slashI = '/';
	}
#endif
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
	if (!Reference->saveToFile(file)) return false;
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
	fprintf(file, "%f %f %u\n", HashTable::MainPara, HashTable::OverflowPara, seed::seedLen);
	fprintf(file, "%u %u %u %u\n", Index->nMainTable, Index->nUsedMainTable, Index->nOverflowTable, Index->nUsedOverflowTable);
	fclose(file);
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
	for (unsigned i = 0; i < Index->nOverflowTable; ++i)
	{
		fwrite(&(Index->OverflowTable[i].key), sizeof(unsigned long long), 1, file);
		unsigned size = Index->OverflowTable[i].locations.size();
		fwrite(&size, sizeof(unsigned), 1, file);
		//fprintf(file, "%llu %u ", Index->OverflowTable[i].key, Index->OverflowTable[i].locations.size());
		if (size == 0) continue;
		fwrite(&(Index->OverflowTable[i].locations[0]), sizeof(std::vector<unsigned>::value_type), size, file);
		//for (unsigned j = 0; j < size; ++j)
		//{
		//	//fprintf(file, "%u ", Index->OverflowTable[i].locations[j]);
		//	fwrite(&(Index->OverflowTable[i].locations[j]), sizeof(std::vector<unsigned>::size_type), 1, file);
		//}
		//fwrite((const void*)&(Index->OverflowTable[i].locations[0]), sizeof(std::vector<unsigned>::size_type), Index->OverflowTable[i].locations.size(), file);
	}
	fclose(file);
	delete tableFileName;
	delete oftFileName;
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
	fscanf(file, "%f %f %u\n", &(HashTable::MainPara), &(HashTable::OverflowPara), &(seed::seedLen));
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
	unsigned nloc;
	Index->OverflowTable = (HashTable::OverflowEntry*)BigAlloc(Index->nOverflowTable*sizeof(HashTable::OverflowEntry));
	for (unsigned i = 0; i < Index->nOverflowTable; ++i)
	{
		//fscanf(file, "%llu %u ", &(Index->OverflowTable[i].key), &nloc);
		fread(&(Index->OverflowTable[i].key), sizeof(unsigned long long), 1, file);
		fread(&nloc, sizeof(unsigned), 1, file);
		//fscanf(file, "%u ", &tmploc);
		if (nloc == 0) continue;
		Index->OverflowTable[i].locations.resize(nloc);
		fread(&(Index->OverflowTable[i].locations[0]), sizeof(std::vector<unsigned>::value_type), nloc, file);
	}
	fclose(file);

	//create the genome filename
	const char *genomeFileExtension = ".genome";
	size_t genomeFileNameLength = (size_t)(strlen(fname) + strlen(genomeFileExtension) + 1);
	char *genomeFileName = new char[genomeFileNameLength];
	_snprintf(genomeFileName, genomeFileNameLength, "%s%s", fname, genomeFileExtension);
	assert(Reference == nullptr);
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
			if (Reference->chrs[i].sequence.isASeed(j, seed::seedLen))
			{
				seed TheSeed = Reference->chrs[i].sequence.getSeed(j, seed::seedLen);
				//seed TheSeed(tmpSeq);
				Index->insert(TheSeed.getBases(), Reference->chrs[i].start_location + j);
			}
		}
	}
	return Index;
}

HashTable* Indexer::cutIndex(unsigned hMax, const char * fname)
{
	if (Index == nullptr || Reference == nullptr) return Index;
	unsigned FinishedIndex = 0;
	if (fname != nullptr)
	{
		saveToFile(fname);
		FILE * file;
		file = fopen(fname, "rb");
		if (file != nullptr)
		{
			fread(&FinishedIndex, sizeof(unsigned), 1, file);
			fclose(file);
		}
	}
	HashTable::Result tmpResult;
	unsigned tmpLocation[2], tmpStartLoc[2], tmpLength[2], tmpChrNum[2], tmpComputeLoc[2], tmpComputeLength[2];
	int tmpEdd;
	unsigned SimilaredCount = 0;
	short(*L)[2 * MAX_K + 1];
	unsigned* PreClear, PreClearCount;
	//#pragma omp parallel private(L,PreClear,PreClearCount)
	{
		L = (short(*)[2 * MAX_K + 1])malloc((MAX_K + 1)*(2 * MAX_K + 1)*sizeof(short));
		PreClear = (unsigned *)malloc((unsigned)(Reference->getTotalLength()*0.01) * sizeof(unsigned));
		for (unsigned block = FinishedIndex; block < Index->nOverflowTable; block += PERSAVE)
		{
			unsigned blockEnd = block + PERSAVE < Index->nOverflowTable ? block + PERSAVE : Index->nOverflowTable;
			//#pragma omp parallel for private(tmpResult, tmpLocation, tmpStartLoc, tmpLength, tmpChrNum, tmpComputeLoc, tmpComputeLength, tmpEdd, SimilaredCount)
			for (long long i = block; i < blockEnd; ++i)
			{
				if (Index->OverflowTable[i].locations.size() > hMax - 2)
				{
					//L = new short[MAX_K + 1][2 * MAX_K + 1];// TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
					SimilaredCount = 0;
					unsigned LocSimilaredCount = 0;
					PreClearCount = 0;
					//Similar location marked as UnusedLocation
					tmpResult = Index->lookup(Index->OverflowTable[i].key);
					for (unsigned x = 0; x < tmpResult.size(); x += 500)
					{
						unsigned xEnd = x + 500 < tmpResult.size() ? x + 500 : tmpResult.size();
						for (unsigned q = x; q < xEnd; ++q)
						{
							LocSimilaredCount = 0;
							tmpLocation[0] = tmpResult[q];
							if (tmpLocation[0] == HashTable::UnusedLocation) continue;
							tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
							tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
							tmpLength[0] = Reference->chrs[tmpChrNum[0]].sequence.length;
							tmpComputeLoc[0] = tmpLocation[0] - ReadLength + seed::seedLen > tmpStartLoc[0] ? tmpLocation[0] - ReadLength + seed::seedLen : tmpStartLoc[0];
							tmpComputeLength[0] = tmpLocation[0] + ReadLength > tmpStartLoc[0] + tmpLength[0] ? tmpStartLoc[0] + tmpLength[0] - tmpComputeLoc[0] : tmpLocation[0] + ReadLength - tmpComputeLoc[0];
							for (unsigned j = q + 1; j < xEnd; ++j)
							{
								tmpLocation[1] = tmpResult[j];
								if (tmpLocation[1] == HashTable::UnusedLocation) continue;
								tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
								tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
								tmpLength[1] = Reference->chrs[tmpChrNum[1]].sequence.length;
								tmpComputeLoc[1] = tmpLocation[1] - ReadLength + seed::seedLen > tmpStartLoc[1] ? tmpLocation[1] - ReadLength + seed::seedLen : tmpStartLoc[1];
								tmpComputeLength[1] = tmpLocation[1] + ReadLength > tmpStartLoc[1] + tmpLength[1] ? tmpStartLoc[1] + tmpLength[1] - tmpComputeLoc[1] : tmpLocation[1] + ReadLength - tmpComputeLoc[1];
								tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.seq + tmpComputeLoc[0] - tmpStartLoc[0], tmpComputeLength[0], Reference->chrs[tmpChrNum[1]].sequence.seq + tmpComputeLoc[1] - tmpStartLoc[1], tmpComputeLength[1], SimilarEdd, L);
								if (tmpEdd >= 0 && tmpEdd <= SimilarEdd)
								{
									tmpResult[j] = HashTable::UnusedLocation;
									++SimilaredCount;
									++LocSimilaredCount;
								}
							}
							if (LocSimilaredCount >= TO_BE_CLEAR)
							{
								PreClear[PreClearCount++] = tmpResult[q];
							}
						}
					}
					//get ride of the unused locations
					double OriginalSize = tmpResult.size();
					if (SimilaredCount != 0)
					{
						unsigned NonVoidCount = 0;
						for (unsigned j = 0; j != tmpResult.size(); ++j)
						{
							if (tmpResult[j] != HashTable::UnusedLocation)
							{
								tmpResult[NonVoidCount++] = tmpResult[j];
							}
						}
						if (SimilaredCount >= Index->OverflowTable[i].locations.size()) Index->OverflowTable[i].locations.clear();
						else  Index->OverflowTable[i].locations.resize(Index->OverflowTable[i].locations.size() - SimilaredCount);
					}
					//estimate if the 2nd cutted size will below hMax, if not we don't do the 2nd cut, because the data will not be used anyway
					if ((unsigned)(tmpResult.size()*(1.0 - (double)SimilaredCount / OriginalSize)) <= hMax)
					{
						SimilaredCount = 0;
						for (unsigned q = 0; q < tmpResult.size(); ++q)
						{
							LocSimilaredCount = 0;
							tmpLocation[0] = tmpResult[q];
							if (tmpLocation[0] == HashTable::UnusedLocation) continue;
							tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
							tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
							tmpLength[0] = Reference->chrs[tmpChrNum[0]].sequence.length;
							tmpComputeLoc[0] = tmpLocation[0] - ReadLength + seed::seedLen > tmpStartLoc[0] ? tmpLocation[0] - ReadLength + seed::seedLen : tmpStartLoc[0];
							tmpComputeLength[0] = tmpLocation[0] + ReadLength > tmpStartLoc[0] + tmpLength[0] ? tmpStartLoc[0] + tmpLength[0] - tmpComputeLoc[0] : tmpLocation[0] + ReadLength - tmpComputeLoc[0];
							for (unsigned j = q + 1; j < tmpResult.size(); ++j)
							{
								tmpLocation[1] = tmpResult[j];
								if (tmpLocation[1] == HashTable::UnusedLocation) continue;
								tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
								tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
								tmpLength[1] = Reference->chrs[tmpChrNum[1]].sequence.length;
								tmpComputeLoc[1] = tmpLocation[1] - ReadLength + seed::seedLen > tmpStartLoc[1] ? tmpLocation[1] - ReadLength + seed::seedLen : tmpStartLoc[1];
								tmpComputeLength[1] = tmpLocation[1] + ReadLength > tmpStartLoc[1] + tmpLength[1] ? tmpStartLoc[1] + tmpLength[1] - tmpComputeLoc[1] : tmpLocation[1] + ReadLength - tmpComputeLoc[1];
								tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.seq + tmpComputeLoc[0] - tmpStartLoc[0], tmpComputeLength[0], Reference->chrs[tmpChrNum[1]].sequence.seq + tmpComputeLoc[1] - tmpStartLoc[1], tmpComputeLength[1], SimilarEdd, L);
								if (tmpEdd >= 0 && tmpEdd <= SimilarEdd)
								{
									tmpResult[j] = HashTable::UnusedLocation;
									++SimilaredCount;
									++LocSimilaredCount;
								}
							}
							if (LocSimilaredCount >= TO_BE_CLEAR)
							{
								PreClear[PreClearCount++] = tmpResult[q];
							}
						}
						//take care of PreClears
						for (unsigned j = 0; j != PreClearCount; ++j)
						{
							for (unsigned c = 0; c != tmpResult.size(); ++c)
							if (tmpResult[c] == PreClear[j])
							{
								if (tmpResult[c] != HashTable::UnusedLocation)//3rd edition change: there might be some location has already be cleared
								{
									tmpResult[c] = HashTable::UnusedLocation;
									++SimilaredCount;
								}
								break;
							}
						}
						//get ride of the cleared ones
						double OriginalSize = tmpResult.size();
						if (SimilaredCount != 0)
						{
							unsigned NonVoidCount = 0;
							for (unsigned j = 0; j != tmpResult.size(); ++j)
							{
								if (tmpResult[j] != HashTable::UnusedLocation)
								{
									tmpResult[NonVoidCount++] = tmpResult[j];
								}
							}
							if (SimilaredCount >= Index->OverflowTable[i].locations.size()) Index->OverflowTable[i].locations.clear();
							else  Index->OverflowTable[i].locations.resize(Index->OverflowTable[i].locations.size() - SimilaredCount);
						}
					}
				}
			}
			std::cout << blockEnd <<" / "<<Index->nOverflowTable<< " overflow entries have been cutted" << std::endl;
			if (fname != nullptr)
			{
				FILE* file;
				file = fopen(fname, "wb");
				if (file == nullptr)
				{
					throw - 102;
				}
				FinishedIndex = blockEnd;
				fwrite(&FinishedIndex, sizeof(unsigned), 1, file);
				fclose(file);
				const char *oftFileExtension = ".overflowtable";
				size_t oftFileNameLength = (size_t)(strlen(fname) + strlen(oftFileExtension) + 1);
				char *oftFileName = new char[oftFileNameLength];
				_snprintf(oftFileName, oftFileNameLength, "%s%s", fname, oftFileExtension);
				//save the Index to the overflow table
				file = fopen(oftFileName, "wb");
				if (file == nullptr)
				{
					throw - 102;
				}
				for (unsigned i = 0; i < Index->nOverflowTable; ++i)
				{
					fwrite(&(Index->OverflowTable[i].key), sizeof(unsigned long long), 1, file);
					unsigned size = Index->OverflowTable[i].locations.size();
					fwrite(&size, sizeof(unsigned), 1, file);
					if (size == 0) continue;
					fwrite(&(Index->OverflowTable[i].locations[0]), sizeof(std::vector<unsigned>::value_type), size, file);
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
						throw - 110;//Ð´ÈëIndexÊ§°Ü£¡
						return false;
					}
					fclose(file);
					++No;
					if (Index->nMainTable <= MAX_SAVE_SIZE*No) break;
				} while (No < 1000);
				delete tableFileName;
			}
		}
		delete L;
		delete PreClear;
	}
	return Index;
}

HashTable* Indexer::throughlyCut(unsigned hMax, const char * fname)
{
	if (Index == nullptr || Reference == nullptr) return Index;
	unsigned FinishedIndex = 0;
	if (fname != nullptr)
	{
		saveToFile(fname);
		FILE * file;
		file = fopen(fname, "rb");
		if (file != nullptr)
		{
			fread(&FinishedIndex, sizeof(unsigned), 1, file);
			fclose(file);
		}
	}
	HashTable::Result tmpResult;
	unsigned tmpLocation[2], tmpStartLoc[2], tmpLength[2], tmpChrNum[2], tmpComputeLoc[2], tmpComputeLength[2];
	int tmpEdd;
	unsigned SimilaredCount = 0;
	short(*L)[2 * MAX_K + 1];
	unsigned* PreClear, PreClearCount;
	//#pragma omp parallel private(L,PreClear,PreClearCount)
	{
		L = (short(*)[2 * MAX_K + 1])malloc((MAX_K + 1)*(2 * MAX_K + 1)*sizeof(short));
		PreClear = (unsigned *)malloc((unsigned)(Reference->getTotalLength()*0.01) * sizeof(unsigned));
		for (unsigned block = FinishedIndex; block < Index->nOverflowTable; block += PERSAVE)
		{
			unsigned blockEnd = block + PERSAVE < Index->nOverflowTable ? block + PERSAVE : Index->nOverflowTable;
			//#pragma omp parallel for private(tmpResult, tmpLocation, tmpStartLoc, tmpLength, tmpChrNum, tmpComputeLoc, tmpComputeLength, tmpEdd, SimilaredCount)
			for (long long i = block; i < blockEnd; ++i)
			{
				if (Index->OverflowTable[i].locations.size() > 0)
				{
					//L = new short[MAX_K + 1][2 * MAX_K + 1];// TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
					SimilaredCount = 0;
					unsigned LocSimilaredCount = 0;
					PreClearCount = 0;
					//Similar location marked as UnusedLocation
					tmpResult = Index->lookup(Index->OverflowTable[i].key);
					for (unsigned x = 0; x < tmpResult.size(); x += 500)
					{
						unsigned xEnd = x + 500 < tmpResult.size() ? x + 500 : tmpResult.size();
						for (unsigned q = x; q < xEnd; ++q)
						{
							LocSimilaredCount = 0;
							tmpLocation[0] = tmpResult[q];
							if (tmpLocation[0] == HashTable::UnusedLocation) continue;
							tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
							tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
							tmpLength[0] = Reference->chrs[tmpChrNum[0]].sequence.length;
							tmpComputeLoc[0] = tmpLocation[0] - ReadLength + seed::seedLen > tmpStartLoc[0] ? tmpLocation[0] - ReadLength + seed::seedLen : tmpStartLoc[0];
							tmpComputeLength[0] = tmpLocation[0] + ReadLength > tmpStartLoc[0] + tmpLength[0] ? tmpStartLoc[0] + tmpLength[0] - tmpComputeLoc[0] : tmpLocation[0] + ReadLength - tmpComputeLoc[0];
							for (unsigned j = q + 1; j < xEnd; ++j)
							{
								tmpLocation[1] = tmpResult[j];
								if (tmpLocation[1] == HashTable::UnusedLocation) continue;
								tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
								tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
								tmpLength[1] = Reference->chrs[tmpChrNum[1]].sequence.length;
								tmpComputeLoc[1] = tmpLocation[1] - ReadLength + seed::seedLen > tmpStartLoc[1] ? tmpLocation[1] - ReadLength + seed::seedLen : tmpStartLoc[1];
								tmpComputeLength[1] = tmpLocation[1] + ReadLength > tmpStartLoc[1] + tmpLength[1] ? tmpStartLoc[1] + tmpLength[1] - tmpComputeLoc[1] : tmpLocation[1] + ReadLength - tmpComputeLoc[1];
								tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.seq + tmpComputeLoc[0] - tmpStartLoc[0], tmpComputeLength[0], Reference->chrs[tmpChrNum[1]].sequence.seq + tmpComputeLoc[1] - tmpStartLoc[1], tmpComputeLength[1], Aligner::confidence, L);
								if (tmpEdd >= 0 && tmpEdd <= Aligner::confidence)
								{
									tmpResult[j] = HashTable::UnusedLocation;
									++SimilaredCount;
									++LocSimilaredCount;
								}
							}
							if (LocSimilaredCount >= TO_BE_CLEAR)
							{
								PreClear[PreClearCount++] = tmpResult[q];
							}
						}
					}
					//get ride of the unused locations
					double OriginalSize = tmpResult.size();
					if (SimilaredCount != 0)
					{
						unsigned NonVoidCount = 0;
						for (unsigned j = 0; j != tmpResult.size(); ++j)
						{
							if (tmpResult[j] != HashTable::UnusedLocation)
							{
								tmpResult[NonVoidCount++] = tmpResult[j];
							}
						}
						if (SimilaredCount >= Index->OverflowTable[i].locations.size()) Index->OverflowTable[i].locations.clear();
						else  Index->OverflowTable[i].locations.resize(Index->OverflowTable[i].locations.size() - SimilaredCount);
					}
					//estimate if the 2nd cutted size will below hMax, if not we don't do the 2nd cut, because the data will not be used anyway
					if ((unsigned)(tmpResult.size()*(1.0 - (double)SimilaredCount / OriginalSize)) <= hMax)
					{
						SimilaredCount = 0;
						for (unsigned q = 0; q < tmpResult.size(); ++q)
						{
							LocSimilaredCount = 0;
							tmpLocation[0] = tmpResult[q];
							if (tmpLocation[0] == HashTable::UnusedLocation) continue;
							tmpChrNum[0] = Reference->getChrNumber(tmpLocation[0]);
							tmpStartLoc[0] = Reference->chrs[tmpChrNum[0]].start_location;
							tmpLength[0] = Reference->chrs[tmpChrNum[0]].sequence.length;
							tmpComputeLoc[0] = tmpLocation[0] - ReadLength + seed::seedLen > tmpStartLoc[0] ? tmpLocation[0] - ReadLength + seed::seedLen : tmpStartLoc[0];
							tmpComputeLength[0] = tmpLocation[0] + ReadLength > tmpStartLoc[0] + tmpLength[0] ? tmpStartLoc[0] + tmpLength[0] - tmpComputeLoc[0] : tmpLocation[0] + ReadLength - tmpComputeLoc[0];
							for (unsigned j = q + 1; j < tmpResult.size(); ++j)
							{
								tmpLocation[1] = tmpResult[j];
								if (tmpLocation[1] == HashTable::UnusedLocation) continue;
								tmpChrNum[1] = Reference->getChrNumber(tmpLocation[1]);
								tmpStartLoc[1] = Reference->chrs[tmpChrNum[1]].start_location;
								tmpLength[1] = Reference->chrs[tmpChrNum[1]].sequence.length;
								tmpComputeLoc[1] = tmpLocation[1] - ReadLength + seed::seedLen > tmpStartLoc[1] ? tmpLocation[1] - ReadLength + seed::seedLen : tmpStartLoc[1];
								tmpComputeLength[1] = tmpLocation[1] + ReadLength > tmpStartLoc[1] + tmpLength[1] ? tmpStartLoc[1] + tmpLength[1] - tmpComputeLoc[1] : tmpLocation[1] + ReadLength - tmpComputeLoc[1];
								tmpEdd = NASeq::computeEditDistance(Reference->chrs[tmpChrNum[0]].sequence.seq + tmpComputeLoc[0] - tmpStartLoc[0], tmpComputeLength[0], Reference->chrs[tmpChrNum[1]].sequence.seq + tmpComputeLoc[1] - tmpStartLoc[1], tmpComputeLength[1], Aligner::confidence, L);
								if (tmpEdd >= 0 && tmpEdd <= Aligner::confidence)
								{
									tmpResult[j] = HashTable::UnusedLocation;
									++SimilaredCount;
									++LocSimilaredCount;
								}
							}
							if (LocSimilaredCount >= TO_BE_CLEAR)
							{
								PreClear[PreClearCount++] = tmpResult[q];
							}
						}
						//take care of PreClears
						for (unsigned j = 0; j != PreClearCount; ++j)
						{
							for (unsigned c = 0; c != tmpResult.size(); ++c)
							if (tmpResult[c] == PreClear[j])
							{
								if (tmpResult[c] != HashTable::UnusedLocation)//3rd edition change: there might be some location has already be cleared
								{
									tmpResult[c] = HashTable::UnusedLocation;
									++SimilaredCount;
								}
								break;
							}
						}
						//get ride of the cleared ones
						double OriginalSize = tmpResult.size();
						if (SimilaredCount != 0)
						{
							unsigned NonVoidCount = 0;
							for (unsigned j = 0; j != tmpResult.size(); ++j)
							{
								if (tmpResult[j] != HashTable::UnusedLocation)
								{
									tmpResult[NonVoidCount++] = tmpResult[j];
								}
							}
							if (SimilaredCount >= Index->OverflowTable[i].locations.size()) Index->OverflowTable[i].locations.clear();
							else  Index->OverflowTable[i].locations.resize(Index->OverflowTable[i].locations.size() - SimilaredCount);
						}
					}
				}
			}
			std::cout << blockEnd << " overflow entries have been cutted" << std::endl;
			if (fname != nullptr)
			{
				FILE* file;
				file = fopen(fname, "wb");
				if (file == nullptr)
				{
					throw - 102;
				}
				FinishedIndex = blockEnd;
				fwrite(&FinishedIndex, sizeof(unsigned), 1, file);
				fclose(file);
				const char *oftFileExtension = ".overflowtable";
				size_t oftFileNameLength = (size_t)(strlen(fname) + strlen(oftFileExtension) + 1);
				char *oftFileName = new char[oftFileNameLength];
				_snprintf(oftFileName, oftFileNameLength, "%s%s", fname, oftFileExtension);
				//save the Index to the overflow table
				file = fopen(oftFileName, "wb");
				if (file == nullptr)
				{
					throw - 102;
				}
				for (unsigned i = 0; i < Index->nOverflowTable; ++i)
				{
					fwrite(&(Index->OverflowTable[i].key), sizeof(unsigned long long), 1, file);
					unsigned size = Index->OverflowTable[i].locations.size();
					fwrite(&size, sizeof(unsigned), 1, file);
					if (size == 0) continue;
					fwrite(&(Index->OverflowTable[i].locations[0]), sizeof(std::vector<unsigned>::value_type), size, file);
				}
				fclose(file);
			}
		}
		delete L;
		delete PreClear;
	}
	return Index;
}