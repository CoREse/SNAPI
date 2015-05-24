#pragma once
#include "Genome.h"
#include "HashTable.h"
#include "defines.h"

//read and save the reference genome, and build and save the hash table
class Indexer
{
public:
	static Genome* Reference;
	static HashTable* Index;
	static unsigned MaxH,ReadLength;
	Genome* readReference(const char * FolderPath);
	HashTable* buildIndex();
	HashTable* cutIndex(unsigned hMax);//the soul of the improved algorithm
	bool saveToFile(const char * FileName);
	bool loadFromFile(const char * FileName);
	void cleanIndex();
	Indexer();
	~Indexer();
};

