#pragma once
#include <stdio.h>
#include <vector>
#include "defines.h"
#include "NASeq.h"
class index
{
public:
	index();
	~index();
	std::vector<NASeq> genomes;
	struct piece
	{
		char * name;
		unsigned offset;
	};
	std::vector<piece> pieces;
	bool readFaFiles(std:: vector<const char*>);
	bool build();//build the hash table
	bool save2file(const char *);
	bool readFromFile(const char *);







	char name[128];//the name of input chromosome(the first line of the input fa file)
	bool build(const char *);//read the sequence from the input file and build the index, return true if succeeded
	bool save2file(const char *);//save to file, return true if succeeded
	bool readFromFile(const char *);//read index, return true if succeeded
	struct entry
	{
		uint64 key;
		unsigned value;
		bool islast;
	} * Index;
	static inline uint64 hash(uint64);
	typedef unsigned size_t;
	static double slack;//the slack parameter
	inline size_t getSize() const { return size; }
	inline size_t getLengthOfGenome() const { return lengthOfGenome; }
	inline const char * getTheGenome() const { return theGenome; }
	size_t lookup(uint64, size_t*&) const;//to lookup the key from the Index, return the number found, and the second parameter return the begining pointer of found entries' indexes, return 0xffffffff means something wrong happened
private:
	size_t size;
	int insert(entry);//insert the entry into Index, return 0:there is already the entry in the Index, 1:successful inserted, -1 wrong
	size_t lengthOfGenome;
	size_t usedIndex;
	char * theGenome;
};

