#pragma once
#include "Genome.h"
#include "HashTable.h"

//read and save the reference genome, and build and save the hash table
class Indexer
{
public:
	static Genome* Reference;
	static HashTable* Index;
	static unsigned MaxH;
	Genome* readReference(const char * FolderPath);
	HashTable* buildIndex();
	bool saveToFile(const char * FileName);
	bool loadFromFile(const char * FileName);
	Indexer();
	~Indexer();
};

