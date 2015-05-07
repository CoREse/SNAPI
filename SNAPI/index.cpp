#include "index.h"
#include "seed.h"
#include <stdlib.h>
#include <string.h>

index::index()
:lengthOfGenome(0), usedIndex(0), Index(nullptr), theGenome(nullptr)
{
}


index::~index()
{
	free(Index);
}

inline uint64 index::hash(uint64 h)
{
	h ^= h >> 33;
	h *= 0xff51afd7ed558ccd;
	h ^= h >> 33;
	h *= 0xc4ceb9fe1a85ec53;
	h ^= h >> 33;
	return h;
}

bool index::build(const char * ifileName)
{
	if (ifileName == nullptr)
	{
		return false;
	}
	FILE *ifile = fopen(ifileName, "r");
	if (ifile == nullptr)
	{
		throw - 102;//file pointer is null
		return false;
	}
	//allocate space
	fseek(ifile, 0, SEEK_END);
	size = (size_t)(double(ftell(ifile))*(1.0 + slack));
	Index = (entry*)malloc(size*sizeof(entry));
	theGenome = (char*)malloc(ftell(ifile)*sizeof(char));
	if (Index == nullptr || theGenome == nullptr)
	{
		fclose(ifile);
		throw - 106;//malloc分配空间失败
		return false;
	}
	rewind(ifile);
	//initializing
	for (size_t i = 0; i != size; ++i)
	{
		Index[i].key = 0;
		Index[i].value = 0xffffffff;
		Index[i].islast = true;
	}
	//read and hash
	char buffer[BUFFER_SIZE];
	size_t read;
	//read the name
	//这里有个很奇怪的问题，当代码为read = (size_t)fread(buffer + seed::length, 1, BUFFER_SIZE, ifile);时，运行时会出错，栈在buffer处有错误什么的。太奇怪了。
	read = (size_t)fread(buffer, 1, BUFFER_SIZE, ifile);
	if (read < 128)
	{
		fclose(ifile);
		throw - 103;
		return false;
	}
	size_t nameLength;
	for (nameLength = 0; nameLength < read; ++nameLength)
	{
		if (buffer[nameLength] == '\n') break;
		if (nameLength >= 127)
		{
			fclose(ifile);
			throw - 104;//name is too long!
			return false;
		}
		name[nameLength] = buffer[nameLength];
	}
	name[nameLength] = '\0';
	if (name[0] == '>')
	{
		for (int i = 1; i <= nameLength; ++i)
		{
			name[i - 1] = name[i];
		}
	}
	--nameLength;
	//read the sequence
	rewind(ifile);
	fread(buffer, 1, nameLength, ifile);
	bool firstRound = true;
	size_t reserved = 2 * seed::length;//每次保留的长度
	while (1)
	{
		read = (size_t)fread(buffer + reserved, 1, BUFFER_SIZE - reserved - 1, ifile);
		if (firstRound)
		{
			if (read < reserved)
			{
				fclose(ifile);
				throw - 103;//the genome reference's length is too short
				return false;
			}
		}
		if (read == 0) break;
		buffer[reserved + read] = '\0';
		for (size_t i = firstRound ? reserved : 0; i < (feof(ifile) ? BUFFER_SIZE : read); ++i)
		{
			if (
				(buffer[i] & 0xDF) != 'A' &&
				(buffer[i] & 0xDF) != 'T' &&
				(buffer[i] & 0xDF) != 'G' &&
				(buffer[i] & 0xDF) != 'C' &&
				(buffer[i] & 0xDF) != 'N'
				)
			{
				continue;
			}
			seed theSeed;
			try
			{
				theSeed.acquireSeq(buffer + i);
			}
			catch (int err)
			{
				if (err == -201)
				{
					int numberofbreak = 0;
					for (int j = 0; j < seed::length - 1; ++j)
					{
						if (buffer[i + j + numberofbreak] == '\n')
						{
							++numberofbreak;
							--j;
							continue;
						}
						theGenome[lengthOfGenome + j] = buffer[i + j + numberofbreak] & 0xDF;//只保留大写的基因字母
					}
					lengthOfGenome += seed::length - 1;
					break;
				}
			}
			theGenome[lengthOfGenome] = theSeed.Seq[0];
			++lengthOfGenome;
			//去除全N的种子（read里不会有全N的种子）
			bool isAllN = true;
			for (int j = 0; j != seed::length; ++j)
			{
				if (theSeed.Seq[j] != 'N')
				{
					isAllN = false;
					break;
				}
			}
			if (isAllN)
			{
				continue;
			}
			//insert the seed into Index
			uint64 key = theSeed.toNum();
			if (key == 3041094791252ui64)
			{
				int a = 0;
			}
			entry ent = { key, lengthOfGenome };//position start at 1
			switch (insert(ent))
			{
			case 1:
				++usedIndex;
				break;
			case -1:
				fclose(ifile);
				throw - 105;//insert fail!
				return false;
			}
		}
		if (feof(ifile)) break;
		for (int i = 0; i != reserved; ++i)
		{
			buffer[i] = buffer[read + i];
		}
		firstRound = false;
	}
	fclose(ifile);
	return true;
}

int index::insert(index::entry ent)
{
	size_t tableIndex = hash(ent.key) % size;
	bool wrapped = false;
	size_t nProbe = 1;
	while (Index[tableIndex].value != 0xffffffff)
	{
		if (Index[tableIndex].key == ent.key)
		{
			if (Index[tableIndex].value == ent.value)
			{
				return 0;
			}
			else
			{
				Index[tableIndex].islast = false;
			}
		}

		if (nProbe <= PROBE_DEPTH)
		{
			tableIndex += nProbe*nProbe;
			++nProbe;
		}
		else
		{
			++tableIndex;
		}
		if (tableIndex >= size)
		{
			if (wrapped)
			{
				return -1;
			}
			tableIndex %= size;
			wrapped = true;
		}
	}
	if (Index[tableIndex].key != 0) return -1;
	Index[tableIndex].key = ent.key;
	Index[tableIndex].value = ent.value;
	return 1;
}

index::size_t index::lookup(uint64 key, index::size_t* &result) const
{
	size_t tableIndex = hash(key) % size;
	size_t nProbe = 1;
	size_t count = 0;
	size_t returnSize = sizeof(size_t);
	result = (size_t*)malloc(returnSize);
	bool wrapped = false;
	bool islast;
	do
	{
		while (Index[tableIndex].key != key&&Index[tableIndex].value != 0xffffffff)
		{
			if (nProbe <= PROBE_DEPTH)
			{
				tableIndex += nProbe*nProbe;
			}
			else
			{
				++tableIndex;
			}
			if (tableIndex >= size)
			{
				if (wrapped)
				{
					break;
				}
				tableIndex %= size;
				wrapped = true;
			}
			++nProbe;
		}
		//no key found.
		if (Index[tableIndex].key != key)
		{
			if (count == 0)
			{
				free(result);
				result = nullptr;
				return count;
			}
			else
			{
				throw - 108;//abnormal islast flag;
				return count;
			}
		}
		if (Index[tableIndex].value == 0xffffffff)
		{
			if (key == 0)
			{
				throw - 108;
				return count;
			}
			throw - 107;//something wrong with the data in the Index
			return 0xffffffff;
		}
		if (count*sizeof(size_t) >= returnSize)
		{
			returnSize += (count > MALLOC_BLOCK ? MALLOC_BLOCK : count)*sizeof(size_t);
			result = (size_t *)realloc(result, returnSize);
			if (result == nullptr)
			{
				free(result);
				throw - 109;//error occured while allocate space
				return 0xffffffff;
			}
		}
		result[count++] = tableIndex;
		islast = Index[tableIndex].islast;
		if (nProbe <= PROBE_DEPTH)
		{
			tableIndex += nProbe*nProbe;
		}
		else
		{
			++tableIndex;
		}
		if (tableIndex >= size)
		{
			if (wrapped)
			{
				break;
			}
			tableIndex %= size;
			wrapped = true;
		}
		++nProbe;
	} while (!islast);
	return count;
}

bool index::save2file(const char * fname)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file = fopen(fname, "w");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	fprintf(file, "%s %u %u %u %u", name, size, lengthOfGenome, usedIndex, seed::length);
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
		if (fwrite(Index + MAX_SAVE_SIZE*No, sizeof(entry), size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE, file) != (size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE))
		{
			fclose(file);
			throw - 110;//写入Index失败！
			return false;
		}
		fclose(file);
		++No;
		if (size <= MAX_SAVE_SIZE*No) break;
	} while (No < 1000);
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
	if (fwrite(theGenome, sizeof(char), lengthOfGenome, file) != lengthOfGenome)
	{
		fclose(file);
		throw - 110;//写入Index失败！
		return false;
	}
	fclose(file);
	return true;
}

bool index::readFromFile(const char * fname)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file = fopen(fname, "r");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	fscanf(file, "%s%u%u%u%u", name, &size, &lengthOfGenome, &usedIndex, &seed::length);
	fclose(file);
	//allocate space
	free(Index);
	Index = (entry*)malloc(size*sizeof(entry));
	if (Index == nullptr)
	{
		fclose(file);
		throw - 106;//malloc分配空间失败
		return false;
	}
	//create the table filenames
	const char *tableFileExtension = ".table";
	size_t tableFileNameLength = (size_t)(strlen(fname) + strlen(tableFileExtension) + 4);
	char *tableFileName = new char[tableFileNameLength];
	int No = 0;
	do
	{
	_snprintf(tableFileName, tableFileNameLength, "%s%s%03d", fname, tableFileExtension, No);
	//read Index from the table files
	file = fopen(tableFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	if (fread(Index + MAX_SAVE_SIZE*No, sizeof(entry), size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No : MAX_SAVE_SIZE, file) != (size <= MAX_SAVE_SIZE*(No + 1) ? size - MAX_SAVE_SIZE*No:MAX_SAVE_SIZE))
	{
		fclose(file);
		throw - 111;//读取Index失败！
		return false;
	}
	fclose(file);
	++No;
	if (size <= MAX_SAVE_SIZE*No) break;
	} while (No < 1000);
	//create the genome filename
	const char *genomeFileExtension = ".genome";
	size_t genomeFileNameLength = (size_t)(strlen(fname) + strlen(genomeFileExtension) + 1);
	char *genomeFileName = new char[genomeFileNameLength];
	_snprintf(genomeFileName, genomeFileNameLength, "%s%s", fname, genomeFileExtension);
	//allocate space
	free(theGenome);
	theGenome = (char*)malloc(lengthOfGenome*sizeof(char));
	if (theGenome == nullptr)
	{
		fclose(file);
		throw - 106;//malloc分配空间失败
		return false;
	}
	//read theGenome from the genome file
	file = fopen(genomeFileName, "rb");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	if (fread(theGenome, sizeof(char), lengthOfGenome, file) != lengthOfGenome)
	{
		fclose(file);
		throw - 111;//读取Index失败！
		return false;
	}
	fclose(file);
	return true;
}

double index::slack = 0.3;