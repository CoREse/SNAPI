#pragma once

#include <stdio.h>
#include "defines.h"
#include "NASeq.h"

class Aligner
{
	void readAndAlignAndWrite(FILE* readFile, FILE* resultFile);//the FILE**s is to ensure every thread uses the same FILE*s
	inline bool writeResult(FILE* ofile, NASeq& TheRead, const char * readName, bool isReversed, unsigned Position, unsigned editDis);
	bool isFirstWrite;
public:
	class LocationHash
	{
		unsigned size;
	public:
		struct Entry
		{
			unsigned BlockLocation;//location/SAME_POS_BLOCK
			unsigned HitsAndIfScored;//the highest bit marked if it is scored(1 means has been), and the lower 31 bits store the number of Hits;
		};
		//static Entry* DefaultValues;//storage the default values, used to quickly renew
		static unsigned Hash(unsigned key)
		{
			// Hash the key.  Use the hash finalizer from the 64 bit MurmurHash3, http://code.google.com/p/smhasher/wiki/MurmurHash3,
			// which is public domain code.
			//

			key ^= key >> 16;
			key *= 0x85ebca6b;
			key ^= key >> 13;
			key *= 0xc2b2ae35;
			key ^= key >> 16;
			return key;
		}
		unsigned getSize() const
		{
			return size;
		}
		void renew();//to prepare the hash table for new use
		Entry* lookup(unsigned BlockLocation);
		Entry * insert(unsigned BlockLocation);
		LocationHash(unsigned Size = 0);
		~LocationHash();
	private:
		Entry* storage;
		unsigned *UsedEntries, UsedEntriesCount;
	};
	unsigned SequenceAligned;
	static unsigned SeedsToTry, dMax, confidence;
	static float SeedsStartPositions[100];//SeedsStartPositions[round] storage the start positions/seedlength of a round
	void singleAlign(const char * readFileName, const char * resultFileName, unsigned thread = 1);
	Aligner();
	~Aligner();
private:
	int alignARead(NASeq & TheRead, unsigned& best, unsigned &secound, unsigned &location,
		unsigned *Candidates, LocationHash * LocationBlocks, short(*L)[63]);
};

