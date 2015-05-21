#include "HashTable.h"
#include <stdlib.h>

#define DEPT_OF_DETECTION 6

double HashTable::OverflowPara = 0.1;
unsigned long long HashTable::UnusedKey = 0xffffffffffffffff;

HashTable::HashTable(unsigned long long NumberOfEntries)
:nMainTable(NumberOfEntries), nOverflowTable(NumberOfEntries*OverflowPara), nUsedMainTable(0), nUsedOverflowTable(0)
{
	if (NumberOfEntries != 0)
	{
		MainTable = (Entry *)malloc(NumberOfEntries*sizeof(Entry));
		OverflowTable = (OverflowEntry *)malloc(nOverflowTable*sizeof(OverflowEntry));
	}
	for (unsigned long long i = 0; i != NumberOfEntries; ++i)
	{
		if (i < nOverflowTable)
		{
			OverflowTable[i].key = UnusedKey;
		}
		MainTable[i].key = UnusedKey;
	}
}

HashTable::~HashTable()
{
	if (MainTable != nullptr)
		delete MainTable;
	if (OverflowTable != nullptr)
		delete OverflowTable;
}

bool HashTable::insert(unsigned long long key, unsigned long long location)
{
	unsigned long long HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(key) + detect*detect)% nMainTable;
		if (MainTable[HashIndex].key == UnusedKey)
		{
			MainTable[HashIndex].key = key;
			MainTable[HashIndex].location = location;
			++nUsedMainTable;
			return true;
		}
		else if (MainTable[HashIndex].key == key)
		{
			detect=0;
			while (detect <= DEPT_OF_DETECTION)
			{
				HashIndex = (Hash(key)+detect*detect) % nOverflowTable;
				if (OverflowTable[HashIndex].key == UnusedKey)
				{
					OverflowTable[HashIndex].key = key;
					OverflowTable[HashIndex].locations.push_back(location);
					++nUsedOverflowTable;
					return true;
				}
				else if (OverflowTable[HashIndex].key == key)
				{
					OverflowTable[HashIndex].locations.push_back(location);
					return true;
				}
				++detect;
			}
		}
		++detect;
	}
	return false;
}
HashTable::Entry* HashTable::lookupMainTable(unsigned long long key)
{
	unsigned long long HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(key) + detect*detect)% nMainTable;
		if (MainTable[HashIndex].key == UnusedKey)
		{
			return nullptr;
		}
		else if (MainTable[HashIndex].key == key)
		{
			return MainTable + HashIndex;
		}
		++detect;
	}
	return nullptr;
}
HashTable::OverflowEntry * HashTable::lookupOverflowTable(unsigned long long key)
{
	unsigned long long HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(key)+detect*detect) % nOverflowTable;
		if (OverflowTable[HashIndex].key == UnusedKey)
		{
			return nullptr;
		}
		else if (OverflowTable[HashIndex].key == key)
		{
			return OverflowTable + HashIndex;
		}
		++detect;
	}
	return nullptr;
}
void HashTable::lookup(unsigned long long key, Entry* MainResult, OverflowEntry* OverflowResult)
{
	MainResult = lookupMainTable(key);
	OverflowResult = lookupOverflowTable(key);
}
