#include "Aligner.h"
#include "NASeq.h"
#include "Indexer.h"
#include "seed.h"
#include <cassert>
#include <thread>
#include <mutex>

#define DEPT_OF_DETECTION 10

char All2s[MAX_SEQUENCE_LENGTH + 1] =
{
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
};

unsigned Aligner::SeedsToTry = 25, Aligner::dMax = 8, Aligner::confidence = 2;
float Aligner::SeedsStartPositions[100];
//Aligner::LocationHash::Entry*  Aligner::LocationHash::DefaultValues = nullptr;

Aligner::LocationHash::LocationHash(unsigned Size)
	:size(Size), UsedEntriesCount(0)
{
	if (size > 0)
	{
		storage = (Entry*)malloc(size*sizeof(Entry));
		UsedEntries = (unsigned*)malloc(size*sizeof(unsigned));
		for (unsigned i = 0; i != size; ++i)
		{
			storage[i].BlockLocation = HashTable::UnusedLocation;
			storage[i].HitsAndIfScored = 0;
		}
	}
}

Aligner::LocationHash::~LocationHash()
{
	free(storage);
	free(UsedEntries);
}

void Aligner::LocationHash::renew()
{
	for (unsigned i = 0; i < UsedEntriesCount; ++i)
	{
		storage[UsedEntries[i]].BlockLocation = HashTable::UnusedLocation;
	}
	UsedEntriesCount = 0;
}

Aligner::LocationHash::Entry* Aligner::LocationHash::lookup(unsigned BlockLocation)
{
	unsigned HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(BlockLocation) + detect*detect) % size;
		if (storage[HashIndex].BlockLocation == HashTable::UnusedLocation)
		{
			return nullptr;
		}
		else if (storage[HashIndex].BlockLocation == BlockLocation)
		{
			return storage + HashIndex;
		}
		++detect;
	}
	return nullptr;
}

Aligner::LocationHash::Entry* Aligner::LocationHash::insert(unsigned BlockLocation)
{
	unsigned HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(BlockLocation) + detect*detect) % size;
		if (storage[HashIndex].BlockLocation == HashTable::UnusedLocation)
		{
			storage[HashIndex].BlockLocation = BlockLocation;
			storage[HashIndex].HitsAndIfScored = 1;
			UsedEntries[UsedEntriesCount++] = HashIndex;
			return storage + HashIndex;
		}
		if (storage[HashIndex].BlockLocation == BlockLocation)
		{
			++(storage[HashIndex].HitsAndIfScored);
		}
		++detect;
	}
	return nullptr;
}

Aligner::Aligner()
	:isFirstWrite(true), SequenceAligned(0)
{
}

Aligner::~Aligner()
{
}

std::vector<std::thread> threads;
std::mutex mtxR, mtxW;

void Aligner::singleAlign(const char * readFileName, const char * resultFileName, unsigned thread)
{
	if (Indexer::Index == nullptr || Indexer::Reference == nullptr)
	{
		return;
	}
	if (readFileName == nullptr)
	{
		return;
	}
	FILE * ReadFile, *SamFile;
	ReadFile = fopen(readFileName, "r");
	if (ReadFile == nullptr)
	{
		throw - 102;
		return;
	}
	if (resultFileName[0] == '\0')
	{
		resultFileName = "result.sam";
	}
	SamFile = fopen(resultFileName, "w");
	if (SamFile == nullptr)
	{
		throw - 102;
		return;
	}
	for (int i = 0; i != thread; ++i)
	{
		threads.push_back(std::thread([&]() {
			readAndAlignAndWrite(ReadFile, SamFile);
		}));
	}
	for (auto& athread : threads) {
		athread.join();
	}
	fclose(ReadFile);
	fclose(SamFile);
}

void Aligner::readAndAlignAndWrite(FILE* readFile, FILE* resultFile)
{
	short L[MAX_K + 1][2 * MAX_K + 1];// TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
	char *readName = new char[BUFFER_SIZE];
	unsigned *Candidates = new unsigned[MAX_SEED_NUM*MAX_HIT_NUM];
	NASeq TheRead, ReversedRead;
	LocationHash LocationBlocks(SeedsToTry*Indexer::hMax * 10);
	unsigned obest, osecond, olocation, rbest, rsecond, rlocation;
	int oreturn, rreturn;
	bool suc = true;
	while (suc)
	{
		mtxR.lock();
		suc = TheRead.readNext(readFile, readName);
		if (!suc)
		{
			mtxR.unlock();
			break;
		}
		++SequenceAligned;
		mtxR.unlock();
		LocationBlocks.renew();
		oreturn = alignARead(TheRead, obest, osecond, olocation, Candidates, &LocationBlocks, L);
		if (oreturn == 0 && obest == 0)
		{
			mtxW.lock();
			writeResult(resultFile, TheRead, readName, false, olocation, obest);
			mtxW.unlock();
			continue;
		}
		ReversedRead = TheRead;
		ReversedRead.reverse();
		rreturn = alignARead(ReversedRead, rbest, rsecond, rlocation, Candidates, &LocationBlocks, L);
		mtxW.lock();
		if (obest < rbest)
		{
			writeResult(resultFile, TheRead, readName, false, olocation, obest);
		}
		else if (rbest < obest)
		{
			writeResult(resultFile, TheRead, readName, true, rlocation, rbest);
		}
		else if (oreturn == 0)
		{
			writeResult(resultFile, TheRead, readName, false, olocation, obest);
		}
		else if (rreturn == 0)
		{
			writeResult(resultFile, TheRead, readName, true, rlocation, rbest);
		}
		else
		{
			writeResult(resultFile, TheRead, readName, false, olocation, obest);
		}
		mtxW.unlock();
	}
	delete readName;
	delete Candidates;
}

bool Aligner::writeResult(FILE* ofile, NASeq& TheRead, const char * readName, bool isReversed, unsigned Position, unsigned editDis)
{
	if (isFirstWrite)
	{
		fprintf(ofile, "@SQ\tSN:");
		for (unsigned i = 0; i != Indexer::Reference->chrs.size(); ++i)
		{
			fprintf(ofile, "%s", Indexer::Reference->chrs[i].name.c_str());
		}
		fprintf(ofile, "\tLN:%u", Indexer::Reference->getTotalLength());
		isFirstWrite = false;
	}
	//char theQuilities[MAX_SEQUENCE_LENGTH + 1];//just some twos
	//for (int i = 0; i != length; ++i)
	//{
	//	theQuilities[i] = '2';
	//}
	//theQuilities[length] = '\0';
	All2s[TheRead.length] = '\0';
	char CIGAR[2] = "-";
	fprintf(ofile, "\n%s\t%d\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t", readName, isReversed ? 16 : 0, Indexer::Reference->chrs[Indexer::Reference->getChrNumber(Position)].name.c_str(), Position + 1, 60, CIGAR, "*", 0, 0);
	fwrite(TheRead.seq, sizeof(char), TheRead.length, ofile);
	fprintf(ofile, "\t%s\tPG:Z:SNAP\tRG:Z:FASTQ", All2s);
	if (Position != 0)
	{
		fprintf(ofile, "\tNM:i:%u", editDis);
	}
	All2s[TheRead.length] = '2';
	return true;
}

//return 0: single hit, 1: multi hits, -1: not found
int Aligner::alignARead(NASeq& TheRead, unsigned& best, unsigned & second, unsigned &location,
	unsigned *Candidates, LocationHash * LocationBlocks, short(*L)[63])
{
	assert(TheRead.length >= seed::seedLen);
	unsigned SeedsCount = 0, Round = 1, LastRoundsLastSeed = 0, FirstRoundSeeds, cn = 0;
	HashTable::Result LocationsResult;
	LocationHash::Entry* tmpEntry;
	best = MAX_EDD;
	second = MAX_EDD;
	while (SeedsCount < SeedsToTry)
	{
		//Get The Seed
		unsigned SeedLoc = seed::seedLen*SeedsStartPositions[Round] + seed::seedLen*(SeedsCount - LastRoundsLastSeed);
		unsigned UnscoredMostHitedLocation = HashTable::UnusedLocation, UnscoredMostHits = 0, UnscoredSecondHitedLocation = HashTable::UnusedLocation, UnscoredSecondMostHits = 0, dLimit;
		if (SeedLoc + seed::seedLen > TheRead.length)
		{
			LastRoundsLastSeed = SeedsCount;
			if (Round == 1) FirstRoundSeeds = SeedsCount;
			SeedLoc = seed::seedLen*SeedsStartPositions[++Round];
		}
		//TheSeedsToTry[SeedsCount] = 
		++SeedsCount;
		if (!TheRead.isASeed(SeedLoc, seed::seedLen))continue;
		seed TheSeed = TheRead.getSeed(SeedLoc, seed::seedLen);
		LocationsResult = Indexer::Index->lookup(TheSeed.getBases());
		//Get the number of Hits
		unsigned nHits = LocationsResult.size();
		if (nHits == 0) continue;
		if (nHits > Indexer::hMax)
		{
#ifdef SNAP
			continue;
#else
			nHits = Indexer::hMax;
#endif
		}
		//clustering the locations and get the most hitted one
		for (unsigned i = 0; i != nHits; ++i)
		{
#ifdef DEBUG
			assert(LocationsResult[i] != HashTable::UnusedLocation);
#endif
			tmpEntry = LocationBlocks->lookup((LocationsResult[i] - SeedLoc) / SAME_POS_BLOCK);
			if (tmpEntry == nullptr)
			{
				Candidates[cn++] = LocationsResult[i] - SeedLoc;
				assert(LocationBlocks->insert((LocationsResult[i] - SeedLoc) / SAME_POS_BLOCK) != nullptr);
				if (UnscoredMostHits < 1)
				{
					UnscoredMostHits = 1;
					UnscoredMostHitedLocation = LocationsResult[i] - SeedLoc;
				}
			}
			else
			{
				++(tmpEntry->HitsAndIfScored);
				if (tmpEntry->HitsAndIfScored < 0x80000000)
				{
					if (UnscoredMostHits < tmpEntry->HitsAndIfScored)
					{
						UnscoredSecondMostHits = UnscoredMostHits;
						UnscoredSecondHitedLocation = UnscoredMostHitedLocation;
						UnscoredMostHits = (tmpEntry->HitsAndIfScored) & 0x7fffffff;
						UnscoredMostHitedLocation = LocationsResult[i] - SeedLoc;
					}
				}
			}
		}
		//score the UnscoredMostHitedLocation
		if (UnscoredMostHitedLocation == HashTable::UnusedLocation) continue;
		if (best > dMax) dLimit = dMax + confidence - 1;
		else if (second >= best + confidence) dLimit = best + confidence - 1;
		else dLimit = best - 1 > 0 ? best - 1 : 0;
		unsigned chrNum = Indexer::Reference->getChrNumber(UnscoredMostHitedLocation);
		int edd = NASeq::computeEditDistance(Indexer::Reference->chrs[chrNum].sequence.seq + UnscoredMostHitedLocation - Indexer::Reference->chrs[chrNum].start_location, TheRead.length, TheRead.seq, TheRead.length, dLimit, L);
		LocationBlocks->lookup(UnscoredMostHitedLocation / SAME_POS_BLOCK)->HitsAndIfScored &= 80000000;
		if (edd > -1)
		{
			if (best > edd)
			{
				second = best;
				best = edd;
				location = UnscoredMostHitedLocation;
			}
			else if (second > edd) second = edd;
		}
		if (best == 0)
		{
			if (second < best + confidence)
				return 1;
		}
		else if (Round == 1 ? (SeedsCount >= best + confidence) : (FirstRoundSeeds >= best + confidence))
		{
			for (unsigned i = 0; i != cn; ++i)
			{
				tmpEntry = LocationBlocks->lookup(Candidates[i] / SAME_POS_BLOCK);
				assert(tmpEntry != nullptr);
				if (tmpEntry->HitsAndIfScored < 0x80000000)
				{
					if (best > dMax) dLimit = dMax + confidence - 1;
					else if (second >= best + confidence) dLimit = best + confidence - 1;
					else dLimit = best - 1 > 0 ? best - 1 : 0;
					chrNum = Indexer::Reference->getChrNumber(Candidates[i]);
					edd = NASeq::computeEditDistance(Indexer::Reference->chrs[chrNum].sequence.seq + Candidates[i] - Indexer::Reference->chrs[chrNum].start_location, TheRead.length, TheRead.seq, TheRead.length, dLimit, L);
					if (edd > -1)
					{
						if (best > edd)
						{
							second = best;
							best = edd;
							location = Candidates[i];
						}
						else if (second > edd) second = edd;
					}
					if (best == 0)
					{
						if (second < best + confidence)
							return 1;
					}
				}
			}
			break;
		}
		UnscoredMostHitedLocation = UnscoredSecondHitedLocation;
		UnscoredMostHits = UnscoredSecondMostHits;
		UnscoredSecondHitedLocation = HashTable::UnusedLocation;
		UnscoredSecondMostHits = 0;
	}
	if (best <= dMax&& second >= best + confidence)
	{
		return 0;
	}
	else if (best <= dMax)
	{
		return 1;
	}
	return -1;
}