#include "NASeq.h"
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "defines.h"

#define MAX_SEQUENCE_LENGTH 3000000000

char buffer[10240];
//char DNAcharset[17] = "VVVVVVVVAGCTNNNN";
char NASeq::DNAcharset[16] =
{
	'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V',//V means vacancy
	'A', 'G', 'C', 'T', 'N', 'N', 'N', 'N'
};
char NASeq::RNAcharset[16] =
{
	'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V',//V means vacancy
	'A', 'G', 'C', 'U', 'N', 'N', 'N', 'N'
};

std::map <base*, unsigned> NASeq::copies;

inline int min(int a, int b)
{
	return a > b ? b : a;
}

short L[MAX_K + 1][2 * MAX_K + 1];// TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells

inline void CountTrailingZeroes(uint64 x, unsigned long &z)
{
	for (z = 1; z <= 64; ++z)
	{
		if (((x >> z) << z) != x)
		{
			--z;
			return;
		}
	}
}

int NASeq::computeEditDistance(const NASeq& a, const NASeq& b, int k)
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
	if (a.length == 0 || b.length == 0) {
		return -1;
	}
	unsigned p = 0, t = 0;
	int end = min(a.length, b.length);
	while (p < end) {
		uint64 xb, xa;
		xb = b.get32bit(p); xa = a.get32bit(t);
		uint64 x = b.get32bit(p) ^ a.get32bit(t);
		if (x) {
			unsigned long zeroes;
			CountTrailingZeroes(x, zeroes);
			zeroes >>= 2;
			L[0][MAX_K] = min((int)(p) + (int)zeroes, end);
			goto done1;
		}
		p += 8;
		t += 8;
	}
	L[0][MAX_K] = end;
done1:
	if (L[0][MAX_K] == end) {
		int result = (b.length >end ? b.length - end : 0); // Could need some deletions at the end, a is the text
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

			unsigned p = best;
			unsigned t = d + best;
			if (b[p].val == a[t].val) {
				int end = min(b.length, a.length - d);

				while (true) {
					uint64 x = b.get32bit(p) ^ a.get32bit(t);
					if (x) {
						unsigned long zeroes;
						CountTrailingZeroes(x, zeroes);
						zeroes >>= 2;
						best = min((int)(p)+(int)zeroes, end);
						break;
					}
					p += 8;
					if (p >= end) {
						best = end;
						break;
					}
					t += 8;
				}
			}

			if (best == b.length) {
				return e;
			}
			L[e][MAX_K + d] = best;
		}
	}
	return -1;
}

NASeq::NASeq()
:isRNA(false), isSkipped(false), LengthInByte(0), length(0), seq(nullptr), originalSeq(nullptr), isSubString(false)
{
}

NASeq::NASeq(const char * sseq, unsigned slength)
: isRNA(false), isSkipped(false), originalSeq(nullptr), isSubString(false)
{
	length = slength == 0 ? strlen(sseq) : slength;
	if (length == 0)
	{
		LengthInByte = 0;
		seq = nullptr;
		return;
	}
	LengthInByte = length / 2 + length % 2;
	seq = (base*)malloc(LengthInByte*sizeof(base));
	for (unsigned i = 0; i < LengthInByte; ++i)
	{
		seq[i] = sseq + i * 2;
	}
	originalSeq = seq;
	copies[originalSeq] = 1;
}

NASeq::NASeq(const NASeq& b)
{
	*this = b;
}

NASeq::NASeq(FILE * ifile, bool isFQ)
:originalSeq(nullptr), isSubString(false)
{
	if (isFQ)
	{
		if (ifile == nullptr)
		{
			throw - 102;
			return;
		}
		fscanf(ifile, "%s", buffer);
		if (feof(ifile)) return;
		if (buffer[0] != '@')
		{
			throw - 112;//file content error
			return;
		}
		fscanf(ifile, "%s", buffer);
		NASeq(buffer);
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
		//allocate space
		length = 0;
		LengthInByte = 0;
		fseek(ifile, 0, SEEK_END);
		//seq = new base[ftell(ifile) / 2 + 1];//the sizeof(base) is 1. The fa file always provide a very long sequence, in this case, we don't mind to allocate a little more space for restoring it
		seq = (base*)malloc(ftell(ifile) / 2 + 1);//oddly, the malloc is far far more faster than the new operater
		if (seq == nullptr)
		{
			fclose(ifile);
			throw - 106;//malloc∑÷≈‰ø’º‰ ß∞‹
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
			if (isSkipped == (length % 2 == 0))
				//ƒ©Œ≤”–ø’
			{
				seq[LengthInByte - 1].val &= 0xf0;
				switch ((*buffer) & 0xDF)
				{
				case 'A':
					seq[LengthInByte - 1].val += 0x08;
					break;
				case 'T':
				case 'U':
					seq[LengthInByte - 1].val += 0x0B;
					break;
				case 'G':
					seq[LengthInByte - 1].val += 0x09;
					break;
				case 'C':
					seq[LengthInByte - 1].val += 0x0a;
					break;
				case 'N':
					seq[LengthInByte - 1].val += 0x0c;
					break;
				default:
					throw 1;
					break;
				}
				for (unsigned i = 0; i != (slen - 1) / 2; ++i)
				{
					seq[LengthInByte + i] = buffer + 1 + i * 2;
				}
				if (slen % 2 == 0)
				{
					seq[LengthInByte + slen / 2 - 1] = buffer[slen - 1];
				}
				LengthInByte += slen / 2;
			}
			else
				//ƒ©Œ≤Œﬁø’
			{
				for (unsigned i = 0; i != slen / 2 + slen % 2; ++i)
				{
					seq[LengthInByte + i] = buffer + 2 * i;
				}
				LengthInByte += slen / 2 + slen % 2;
			}
			length += slen;
		}
	}
	originalSeq = seq;
	copies[originalSeq] = 1;
}

NASeq::~NASeq()
{
	if (!isSubString)
	if (--copies[originalSeq] == 0)
	{
		delete (base*)originalSeq;
		copies.erase(originalSeq);
	}
}

NASeq& NASeq::operator=(const NASeq & b)
{
	if (b.isSubString)
	{
		length = b.length;
		LengthInByte = b.LengthInByte;
		seq = b.seq;
		isSkipped = b.isSkipped;
		isSubString = b.isSubString;
		return (*this);
	}
	++copies[b.originalSeq];
	originalSeq = b.originalSeq;
	isRNA = b.isRNA;
	isSkipped = b.isSkipped;
	length = b.length;
	LengthInByte = b.LengthInByte;
	seq = b.seq;
	isSubString = b.isSubString;
	return *this;
}

/*the old version
NASeq& NASeq::operator=(NASeq & b)
{
++b.copies;
isRNA = b.isRNA;
isSkipped = b.isSkipped;
length = b.length;
LengthInByte = b.LengthInByte;
seq = (base*)malloc(LengthInByte*sizeof(base));
for (unsigned i = 0; i != LengthInByte; ++i)
{
seq[i] = b.seq[i];
}
copies = b.copies;
return *this;
}
*/

NASeq NASeq::operator~()
{
	return reverse();
}

NASeq& NASeq::reverse()
{
	for (unsigned i = 0; i != LengthInByte; ++i)
	{
		~seq[i];
	}
	return *this;
}

/*
NASeq& NASeq::append(const NASeq& b)
{
if (b.length == 0)
return *this;
if (isSkipped == (length % 2 == 0))
//ƒ©Œ≤”–ø’
{
seq = (base*)realloc(seq, LengthInByte + b.length / 2);//sizeof(base) is 1
if (b.isSkipped)
{
seq[LengthInByte - 1].val &= 0xf0;
seq[LengthInByte - 1].val |= b.seq[0].val & 0x0f;
for (unsigned i = 1; i < b.LengthInByte; ++i)
{
seq[LengthInByte + i - 1] = b.seq[i];
}
}
else
{
seq[LengthInByte - 1].val &= 0xf0;
seq[LengthInByte - 1].val |= b.seq[0].val >> 4;
for (unsigned i = 0; i != b.LengthInByte - 1; ++i)
{
seq[LengthInByte + i].val = b.seq[i].val << 4 + b.seq[i + 1].val >> 4;
}
if (b.length % 2 == 0)
{
seq[LengthInByte + b.LengthInByte - 1].val = b.seq[b.LengthInByte - 1].val << 4;
}
}
LengthInByte += b.length / 2;
}
else
//ƒ©Œ≤Œﬁø’
{
seq = (base*)realloc(seq, LengthInByte + b.length / 2 + b.length % 2);//sizeof(base) is 1
if (b.isSkipped)
{
for (unsigned i = 0; i != b.LengthInByte - 1; ++i)
{
seq[LengthInByte + i].val = b.seq[i].val << 4 | (b.seq[i + 1].val >> 4);
}
if (b.length % 2 == 1)
{
seq[LengthInByte + b.LengthInByte - 1].val = b.seq[LengthInByte - 1].val << 4;
}
}
else
{
for (unsigned i = 0; i != b.LengthInByte; ++i)
{
seq[LengthInByte + i] = b.seq[i];
}
}
LengthInByte += b.length / 2 + b.length % 2;
}
length += b.length;
return *this;
}

NASeq& NASeq::append(const char * sseq, unsigned slen)
{
if (slen == 0)
{
slen = strlen(sseq);
}
if (slen == 0)
return *this;
if (isSkipped == (length % 2 == 0))
//ƒ©Œ≤”–ø’
{
seq = (base*)realloc(seq, LengthInByte + slen / 2);//sizeof(base) is 1
seq[LengthInByte - 1].val &= 0xf0;
switch ((*sseq) & 0xDF)
{
case 'A':
seq[LengthInByte - 1].val += 0x08;
break;
case 'T':
case 'U':
seq[LengthInByte - 1].val += 0x0B;
break;
case 'G':
seq[LengthInByte - 1].val += 0x09;
break;
case 'C':
seq[LengthInByte - 1].val += 0x0a;
break;
case 'N':
seq[LengthInByte - 1].val += 0x0c;
break;
default:
throw 1;
break;
}
for (unsigned i = 0; i != (slen - 1) / 2; ++i)
{
seq[LengthInByte + i] = sseq + 1 + i * 2;
}
if (slen % 2 == 0)
{
seq[LengthInByte + slen / 2 - 1] = sseq[slen - 1];
}
LengthInByte += slen / 2;
}
else
//ƒ©Œ≤Œﬁø’
{
seq = (base*)realloc(seq, LengthInByte + slen / 2 + slen % 2);//sizeof(base) is 1
for (unsigned i = 0; i != slen / 2 + slen % 2; ++i)
{
seq[LengthInByte + i] = sseq + 2 * i;
}
LengthInByte += slen / 2 + slen % 2;
}
length += slen;
return *this;
}
*/

/*
char * NASeq::toString(unsigned pos, unsigned length)
{
if (length == 0)
{
length = this->length - pos;
}
//#################################################
//note: need futher job to handle the memory
//####################################################
char * str = new char[length + 1];
if (isRNA)
{
if (isSkipped)
{
for (unsigned i = 0; i != length; ++i)
{
if (pos + i >= this->length) break;
str[i] = RNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val & 0x0F : seq[(pos + i) / 2 + 1].val >> 4];
}
}
else
{
for (unsigned i = 0; i != length; ++i)
{
if (pos + i >= this->length) break;
str[i] = RNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val >> 4 : seq[(pos + i) / 2].val & 0x0F];
}
}
}
else
{
if (isSkipped)
{
for (unsigned i = 0; i != length; ++i)
{
if (pos + i >= this->length) break;
str[i] = DNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val & 0x0F : seq[(pos + i) / 2 + 1].val >> 4];
}
}
else
{
for (unsigned i = 0; i != length; ++i)
{
if (pos + i >= this->length) break;
str[i] = DNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val >> 4 : seq[(pos + i) / 2].val & 0x0F];
}
}
}
str[length] = '\0';
return str;
}
*/

base NASeq::operator[] (unsigned i) const
{
	if (isSkipped)
	{
		return base(seq[i / 2 + i % 2].val << (i % 2 == 0 ? 4 : 0));
	}
	else
	{
		return base(seq[i / 2].val << (i % 2 == 0 ? 0 : 4));
	}
}


NASeq NASeq::getSubSequence(int start, unsigned length)
{
	assert(length != 0 && this->length != 0);
	if (start<0) start = 0;
	unsigned start_loc = start;
	NASeq sub;
	//sub.isRNA = isRNA;
	sub.length = this->length - start_loc > length ? length : this->length - start_loc;
	//sub.originalSeq = originalSeq;
	//++copies[originalSeq];
	sub.isSubString = true;
	if (isSkipped)
	{
		sub.seq = seq + sub.length / 2 + sub.length % 2;
		sub.isSkipped = start_loc % 2 == 0;
		sub.LengthInByte = (sub.isSkipped ? 1 : 0) + sub.length / 2 + (sub.isSkipped ? 0 : sub.length % 2);
	}
	else
	{
		sub.seq = seq + start_loc / 2;
		sub.isSkipped = start_loc % 2;
		sub.LengthInByte = (sub.isSkipped ? 1 : 0) + sub.length / 2 + (sub.isSkipped ? 0 : sub.length % 2);
	}
	return sub;
}


bool NASeq::saveToFile(FILE* file)
{
	fprintf(file, "%d %d %u %u\n", isSkipped, isRNA, length, LengthInByte);
	if (fwrite(seq, sizeof(base), LengthInByte, file) != LengthInByte)
	{
		fclose(file);
		throw - 110;//–¥»ÎNASeq ß∞‹£°
		return false;
	}
	return true;
}

bool NASeq::loadFromFile(FILE* file)
{
	if (seq != nullptr) return false;
	fscanf(file, "%d %d %u %u\n", &isSkipped, &isRNA, &length, &LengthInByte);
	seq = (base*)malloc(LengthInByte);
	originalSeq = seq;
	copies[originalSeq] = 1;
	isSubString = false;
	if (fread(seq, sizeof(base), LengthInByte, file) != LengthInByte)
	{
		fclose(file);
		throw - 110;//–¥»ÎNASeq ß∞‹£°
		return false;
	}
	return true;
}

unsigned long long NASeq::get64bit(unsigned index)const
{
	if (isSkipped)
	{
		if (index % 2 == 0)
		{
			unsigned long long tmp = *((unsigned long long *)(seq + index / 2));
			tmp << 4;
			if (index + 16 > length) return tmp;
			else
			{
				tmp |= ((seq + (index / 2 + 8))->val) >> 4;
				return tmp;
			}
		}
		else
		{
			return *((unsigned long long *)(seq + index / 2 + 1));
		}
	}
	else
	{
		if (index % 2 == 0)
		{
			return *((unsigned long long *)(seq + index / 2));
		}
		else
		{
			unsigned long long tmp = *((unsigned long long *)(seq + index / 2));
			tmp << 4;
			if (index + 16 > length) return tmp;
			else
			{
				tmp |= ((seq + (index / 2 + 8))->val) >> 4;
				return tmp;
			}
		}
	}
}

unsigned long long NASeq::get32bit(unsigned index)const
{
	if (isSkipped)
	{
		if (index % 2 == 0)
		{
			unsigned long long tmp = *((unsigned*)(seq + index / 2));
			tmp << 4;
			if (index + 8 > length) return tmp;
			else
			{
				tmp |= ((seq + (index / 2 + 4))->val) >> 4;
				return tmp;
			}
		}
		else
		{
			return *((unsigned long long *)(seq + index / 2 + 1));
		}
	}
	else
	{
		if (index % 2 == 0)
		{
			return *((unsigned long long *)(seq + index / 2));
		}
		else
		{
			unsigned long long tmp = *((unsigned long long *)(seq + index / 2));
			tmp << 4;
			if (index + 8 > length) return tmp;
			else
			{
				tmp |= ((seq + (index / 2 + 4))->val) >> 4;
				return tmp;
			}
		}
	}
}