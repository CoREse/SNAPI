#pragma once
#include "defines.h"
#include "index.h"
#include <stdio.h>
//You MUST fclose the outputFile yourself when you finished the call if you've writen the sam file
class sequence
{
public:
	sequence();
	~sequence();
	static size_t SeedsToTry;
	static size_t dMax;
	static size_t c;//confidence threshold
	static size_t hMax;
	static bool doWrite;//to determin whether write into sam file or not
	static size_t SequenceAligned;
	static char samName[256];
	typedef unsigned size_t;
	bool readNext(FILE*);
	void rc(sequence&);//acquire the reversed complemented sequence
	bool write(FILE*, index&);
	bool align(index&, sequence*&);
	void writeResult(index&);
	static int computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k);
	static FILE* outputFile;//for close
private:
	char Seq[MAX_SEQUENCE_LENGTH + 1];//it is a c-string, need to be all upper-cased
	char name[256];//the fq read name
	bool isReversed;
	size_t editDis;
	size_t Position;
	size_t length;
	char CIGAR[256];
	FILE* writeHead(const char *, index&);
	int singleAlign(index&, size_t &best, size_t &second, size_t &pos);//return value:-1 fail;0:confident;1: ambiguously
};

