#include "HashTable.h"
#include <stdlib.h>
#include <cassert>
#include "BigAlloc.h"

#define DEPT_OF_DETECTION 6

float HashTable::OverflowPara = 0.015, HashTable::MainPara=0.85;
unsigned long long HashTable::UnusedKey = 0xffffffffffffffff;
unsigned HashTable::UnusedLocation = 0xffffffff;

unsigned moreThanOneLocation;
HashTable::HashTable(unsigned NumberOfEntries)
:nMainTable(NumberOfEntries*MainPara), nOverflowTable(NumberOfEntries*OverflowPara), nUsedMainTable(0), nUsedOverflowTable(0)
{
	if (NumberOfEntries != 0)
	{
		unsigned long long tmp = (unsigned long long)NumberOfEntries*sizeof(Entry);
		MainTable = (Entry *)BigAlloc((unsigned long long)nMainTable*sizeof(Entry));
		OverflowTable = (OverflowEntry *)BigAlloc((unsigned long long)nOverflowTable*sizeof(OverflowEntry));
	}
	for (unsigned i = 0; i != nMainTable; ++i)
	{
		if (i < nOverflowTable)
		{
			OverflowTable[i].key = UnusedKey;
		}
		MainTable[i].key = UnusedKey;
		MainTable[i].location1 = UnusedLocation;
		MainTable[i].location2 = UnusedLocation;
	}
}

HashTable::~HashTable()
{
	if (MainTable != nullptr)
		delete MainTable;
	if (OverflowTable != nullptr)
		delete OverflowTable;
}

bool HashTable::insert(unsigned long long key, unsigned location)
{
	unsigned long long HashIndex;
	unsigned detect = 0;
	while (detect <= DEPT_OF_DETECTION)
	{
		HashIndex = (Hash(key) + detect*detect)% nMainTable;
		if (MainTable[HashIndex].key == UnusedKey)
		{
			MainTable[HashIndex].key = key;
			MainTable[HashIndex].location1 = location;
			++nUsedMainTable;
			return true;
		}
		else if (MainTable[HashIndex].key == key)
		{
			if (MainTable[HashIndex].location2 == UnusedLocation)
			{
				MainTable[HashIndex].location2 = location;
				++moreThanOneLocation;
				return true;
			}
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

bool HashTable::saveToFile(FILE* file)
{
	fprintf(file, "%f %f\n", MainPara, OverflowPara);
	fprintf(file,"%u %u %u %u\n",nMainTable,nUsedMainTable,nOverflowTable,nUsedOverflowTable);
	if (fwrite(MainTable,sizeof(Entry),nMainTable,file)!=nMainTable)
	{
		fclose(file);
		throw -110;
	}
	for (unsigned i=0;i<nOverflowTable;++i)
	{
		fprintf(file,"%llu %u ",OverflowTable[i].key,OverflowTable[i].locations.size());
		for (unsigned j=0;j<OverflowTable[i].locations.size();++j)
		{
			fprintf(file,"%u ",OverflowTable[i].locations[j]);
		}
	}
}

bool HashTable::loadFromFile(FILE* file)
{
	assert(MainTable == nullptr&&OverflowTable == nullptr);
	fscanf(file, "%f %f\n%u %u %u %u\n", &MainPara, &OverflowPara, &nMainTable, &nUsedMainTable, &nOverflowTable, &nUsedOverflowTable);
	MainTable = (Entry*)BigAlloc((unsigned long long)nMainTable*sizeof(Entry));
	if (fread(MainTable,sizeof(Entry),nMainTable,file)!=nMainTable)
	{
		fclose(file);
		throw -110;
	}
	unsigned nloc,tmploc;
	OverflowTable = (OverflowEntry*)BigAlloc((unsigned long long)nOverflowTable*sizeof(OverflowEntry));
	for (unsigned i=0;i<nOverflowTable;++i)
	{
		fscanf(file,"%llu %u ",&(OverflowTable[i].key),&nloc);
		for (unsigned j=0;j<nloc;++j)
		{
			fscanf(file,"%u ",&tmploc);
			OverflowTable[i].locations.push_back(tmploc);
		}
	}
}
