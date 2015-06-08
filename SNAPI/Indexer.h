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
	static unsigned hMax, ReadLength, SimilarEdd;
	Genome* readReference(const char * FolderPath);
	HashTable* buildIndex();
	HashTable* cutIndex(unsigned hMax, const char * FileName = nullptr);//the soul of the improved algorithm
	HashTable* throughlyCut(unsigned hMax, const char * FileName = nullptr);
	bool saveToFile(const char * FileName);
	bool loadFromFile(const char * FileName);
	void cleanIndex();
	Indexer();
	~Indexer();
};

