#include "NASeq.h"
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <intrin.h>
#include "defines.h"
#include "seed.h"
#include <thread>

#define MAX_SEQUENCE_LENGTH 3000000000
#define CountTrailingZeroes(x, ans) {_BitScanForward64(&ans, x);}
#define min(a,b)            (((a) < (b)) ? (a) : (b))

char NASeq::buffer[10240];
//char DNAcharset[17] = "VVVVVVVVAGCTNNNN";

int NASeq::LegibleSeed[256];
//=
//{
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//};

char NASeq::Opposite[256];
//=
//{
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	'T', ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	't', ' ', 'g', ' ', ' ', ' ', 'c', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'a', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
//	' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '
//};

unsigned NASeq::NANumber[256];
//=
//{
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//};

std::map <char*, unsigned> NASeq::Copies;

//inline int min(int a, int b)
//{
//	return a > b ? b : a;
//}

//inline void CountTrailingZeroes(uint64 x, unsigned long &z)
//{
//	for (z = 1; z <= 64; ++z)
//	{
//		if (((x >> z) << z) != x)
//		{
//			--z;
//			return;
//		}
//	}
//}

int NASeq::computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k, short L[MAX_K + 1][2 * MAX_K + 1])
{
	for (int i = 0; i < MAX_K + 1; i++) {
		for (int j = 0; j < 2 * MAX_K + 1; j++) {
			L[i][j] = -2;
		}
	}
	/*if (k >= MAX_K)
	{
	throw - 113;
	return MAX_K;
	}*/
	k = min(MAX_K - 1, k); // enforce limit even in non-debug builds
	if (NULL == text) {
		// This happens when we're trying to read past the end of the genome.
		return -1;
	}
	const char* p = pattern;
	const char* t = text;
	int end = min(patternLen, textLen);
	const char* pend = pattern + end;
	while (p < pend) {
		uint64 x = *((const uint64*)p) ^ *((const uint64*)t);
		if (x) {
			unsigned long zeroes;
			CountTrailingZeroes(x, zeroes);
			zeroes >>= 3;
			L[0][MAX_K] = min((int)(p - pattern) + (int)zeroes, end);
			goto done1;
		}
		p += 8;
		t += 8;
	}
	L[0][MAX_K] = end;
done1:
	if (L[0][MAX_K] == end) {
		int result = (patternLen > end ? patternLen - end : 0); // Could need some deletions at the end
		return result;
	}

	for (int e = 1; e <= k; e++) {
		// Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
		for (int d = 0; d != e + 1; d = (d > 0 ? -d : -d + 1)) {
			int best = L[e - 1][MAX_K + d] + 1; // up
			int left = L[e - 1][MAX_K + d - 1];
			if (left > best)
				best = left;
			int right = L[e - 1][MAX_K + d + 1] + 1;
			if (right > best)
				best = right;

			const char* p = pattern + best;
			const char* t = (text + d) + best;
			if (*p == *t) {
				int end = min(patternLen, textLen - d);
				const char* pend = pattern + end;

				while (true) {
					uint64 x = *((const uint64*)p) ^ *((const uint64*)t);
					if (x) {
						unsigned long zeroes;
						CountTrailingZeroes(x, zeroes);
						zeroes >>= 3;
						best = min((int)(p - pattern) + (int)zeroes, end);
						break;
					}
					p += 8;
					if (p >= pend) {
						best = end;
						break;
					}
					t += 8;
				}
			}

			if (best == patternLen) {
				return e;
			}
			L[e][MAX_K + d] = best;
		}
	}
	return -1;
}

NASeq::NASeq()
:length(0), seq(nullptr), isSubString(false), isQuickCopyable(false), OriginalSeq(nullptr)
{
}

NASeq::NASeq(const char * sseq, unsigned slength, bool isQCB)
: isSubString(false), OriginalSeq(nullptr), isQuickCopyable(isQCB)
{
	length = slength == 0 ? strlen(sseq) : slength;
	if (length == 0) return;
	seq = (char*)malloc(length*sizeof(char));
	memcpy(seq, sseq, length*sizeof(char));
	unsigned i;
	if (length > 8)
	for (i = 0; i < length - 8; i += 8)
	{
		*((unsigned long long*)(seq + i)) &= 0xdfdfdfdfdfdfdfdf;
	}
	for (; i < length; ++i)
	{
		seq[i] &= 0xdf;
	}

	OriginalSeq = seq;
	if (isQuickCopyable)
	{
		Copies[OriginalSeq] = 1;
	}
}

NASeq::NASeq(const NASeq& b, bool isQCB)
{
	*this = b;
	isQuickCopyable = isQCB;
}

NASeq::NASeq(FILE * ifile, bool isFQ, char * seqName)
: isSubString(false), seq(nullptr), OriginalSeq(nullptr)
{
	if (isFQ)
	{
		if (ifile == nullptr)
		{
			throw - 102;
			return;
		}
		isQuickCopyable = false;
		fscanf(ifile, "%s", buffer);
		if (feof(ifile)) return;
		if (buffer[0] != '@')
		{
			throw - 112;//file content error
			return;
		}
		if (seqName != nullptr)
		{
			buffer[MAX_NAME_LENGTH - 1] = '\0';
			strcpy(seqName, buffer + 1);
		}
		fscanf(ifile, "%s", buffer);
		length = strlen(buffer);
		if (length == 0)
		{
			seq = nullptr;
			return;
		}
		seq = (char*)malloc(length*sizeof(char));
		memcpy(seq, buffer, length*sizeof(char));
		unsigned i;
		if (length > 8)
		for (i = 0; i < length - 8; i += 8)
		{
			*((unsigned long long*)(seq + i)) &= 0xdfdfdfdfdfdfdfdf;
		}
		for (; i < length; ++i)
		{
			seq[i] &= 0xdf;
		}

		fscanf(ifile, "%s", buffer);
		fscanf(ifile, "%s", buffer);
	}
	else
	{
		if (ifile == nullptr)
		{
			throw - 102;//file pointer is null
			return;
		}
		isQuickCopyable = true;
		//allocate space
		length = 0;
		fseek(ifile, 0, SEEK_END);
		//seq = new base[ftell(ifile) / 2 + 1];//the sizeof(base) is 1. The fa file always provide a very long sequence, in this case, we don't mind to allocate a little more space for restoring it
		seq = (char*)malloc(ftell(ifile) + 1);//oddly, the malloc is far far more faster than the new operater
		if (seq == nullptr)
		{
			fclose(ifile);
			throw - 106;//malloc·ÖÅä¿Õ¼äÊ§°Ü
			return;
		}
		rewind(ifile);
		//read the sequence
		fscanf(ifile, "%s", buffer);//the first line is name of the sequence, which in this class we don't care
		while (true)
		{
			fscanf(ifile, "%s", buffer);
			if (feof(ifile)) return;
			//append(buffer);//very slow
			unsigned slen = strlen(buffer);
			if (slen == 0)
				continue;
			memcpy(seq + length, buffer, slen*sizeof(char));
			unsigned i;
			if (slen > 8)
			for (i = 0; i < slen - 8; i += 8)
			{
				*((unsigned long long*)(seq + length + i)) &= 0xdfdfdfdfdfdfdfdf;
			}
			for (; i < slen; ++i)
			{
				seq[i + length] &= 0xdf;
			}

			length += slen;
		}
	}

	OriginalSeq = seq;
	if (isQuickCopyable)
	{
		Copies[OriginalSeq] = 1;
	}
}

NASeq::~NASeq()
{
	if (seq != nullptr)
	{
		if (isQuickCopyable)
		{
			if (!isSubString)
			if (--Copies[OriginalSeq] == 0)
			{
				delete OriginalSeq;
				Copies.erase(OriginalSeq);
			}
		}
		else
		{
			delete seq;
		}
	}
}

NASeq& NASeq::operator=(const NASeq & b)
{
	if (b.seq == nullptr)
	{
		seq = nullptr;
		return *this;
	}
	if (b.isSubString)
	{
		length = b.length;
		seq = b.seq;
		isSubString = b.isSubString;
		return (*this);
	}
	if (b.isQuickCopyable)
	{
		++Copies[b.OriginalSeq];
		OriginalSeq = b.OriginalSeq;
		length = b.length;
		seq = b.seq;
		isSubString = b.isSubString;
		isQuickCopyable = b.isQuickCopyable;
	}
	else
	{
		if (length >= b.length&&!isQuickCopyable)
		{
			length = b.length;
			memcpy(seq, b.seq, length);
			return *this;
		}
		length = b.length;
		seq = (char*)malloc(length);
		memcpy(seq, b.seq, length);
		OriginalSeq = seq;
		isSubString = b.isSubString;
		isQuickCopyable = b.isQuickCopyable;
	}
	return *this;
}

NASeq NASeq::operator~()
{
	if (isQuickCopyable)
	{
		isQuickCopyable = false;
		NASeq rOne(*this);
		isQuickCopyable = true;
		rOne.reverse();
		return rOne;
	}
	else
	{
		NASeq rOne(*this);
		rOne.reverse();
		return rOne;
	}
}

NASeq& NASeq::reverse()
{
	char tmp;
	for (unsigned i = 0; i != length / 2; ++i)
	{
		tmp = seq[i];
		seq[i] = Opposite[seq[length - i - 1]];
		seq[length - i - 1] = Opposite[tmp];
	}
	if (length % 2 != 0)
	{
		seq[length / 2 + 1] = Opposite[seq[length / 2 + 1]];
	}
	return *this;
}

NASeq NASeq::getSubSequence(int start, unsigned length)
{
	assert(length != 0 && this->length != 0);
	if (!isQuickCopyable)
	{
		return NASeq();
	}
	if (start<0) start = 0;
	unsigned start_loc = start;
	NASeq sub;
	//sub.isRNA = isRNA;
	sub.length = this->length - start_loc > length ? length : this->length - start_loc;
	//sub.originalSeq = originalSeq;
	//++copies[originalSeq];
	sub.isSubString = true;

	if (isQuickCopyable)
	{
		sub.seq = seq;
	}

	return sub;
}

bool NASeq::saveToFile(FILE* file)
{
	fprintf(file, "%d %u\n", isQuickCopyable, length);
	if (fwrite(seq, sizeof(char), length, file) != length)
	{
		fclose(file);
		throw - 110;//Ð´ÈëNASeqÊ§°Ü£¡
		return false;
	}
	return true;
}

bool NASeq::loadFromFile(FILE* file)
{
	if (seq != nullptr) return false;
	fscanf(file, "%d %u\n", &isQuickCopyable, &length);
	seq = (char*)malloc(length);
	OriginalSeq = seq;
	if (isQuickCopyable)
		Copies[OriginalSeq] = 1;
	isSubString = false;
	if (fread(seq, sizeof(char), length, file) != length)
	{
		fclose(file);
		throw - 110;//Ð´ÈëNASeqÊ§°Ü£¡
		return false;
	}
	return true;
}

unsigned long long NASeq::getSeed(unsigned start, unsigned SeedLength) const
{
	if (start + SeedLength > length) throw - 201;
	if (SeedLength == 0) SeedLength = length - start;
	if (SeedLength > 32) SeedLength = 32;
	unsigned long long result = 0;
	for (unsigned i = 0; i != SeedLength; ++i)
	{
		result <<= 2;
		result += NANumber[seq[start + i]];
	}
	return result;
}

bool NASeq::isASeed(unsigned start, unsigned SeedLength) const
{
	if (start + SeedLength > length)
	{
		return false;
	}
	if (SeedLength == 0)
	{
		SeedLength = length - start;
	}
	if (SeedLength > 32) SeedLength = 32;
	for (int i = 0; i != SeedLength; ++i)
	{
		if (LegibleSeed[seq[i + start]] == 0) return false;
	}
	return true;
}

bool NASeq::readNext(FILE* ifile, char * seqName)
{
	if (ifile == nullptr)
	{
		throw - 102;
		return false;
	}
	isQuickCopyable = false;
	fscanf(ifile, "%s", buffer);
	if (feof(ifile)) return false;
	if (buffer[0] != '@')
	{
		throw - 112;//file content error
		return false;
	}
	if (seqName != nullptr)
	{
		buffer[MAX_NAME_LENGTH - 1] = '\0';
		strcpy(seqName, buffer + 1);//the first character is '@'
	}
	fscanf(ifile, "%s", buffer);
	unsigned bufflen = strlen(buffer);
	if (length >= bufflen)
	{
		memcpy(seq, buffer, bufflen*sizeof(char));
		fscanf(ifile, "%s", buffer);
		fscanf(ifile, "%s", buffer);
		return true;
	}
	length = bufflen;
	if (length == 0)
	{
		seq = nullptr;
		return true;
	}
	seq = (char*)malloc(length*sizeof(char));
	memcpy(seq, buffer, length*sizeof(char));

	fscanf(ifile, "%s", buffer);
	fscanf(ifile, "%s", buffer);
	return true;
}