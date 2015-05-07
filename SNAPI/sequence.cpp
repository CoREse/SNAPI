#include "sequence.h"
#include "RSeed.h"
#include <string.h>
#include <stdlib.h>

char All2s[MAX_SEQUENCE_LENGTH + 1] =
{
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
	"22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222"
};

sequence::sequence()
:isReversed(false), length(0), Position(0)
{
	CIGAR[0] = '-';
	CIGAR[1] = '\0';
}


sequence::~sequence()
{
}

void sequence::rc(sequence& b)
{
	b.isReversed = !isReversed;
	b.length = length;
	int count;
	for (count = 0; count <= MAX_SEQUENCE_LENGTH; ++count)
	{
		if (Seq[count] == '\0') break;
	}
	for (int i = 0; i != count; ++i)
	{
		if (Seq[count - 1 - i] == 'N')
		{
			b.Seq[i] = 'N';
		}
		else
		{
			b.Seq[i] = bases[seed::convertBase2Num(Seq[count - 1 - i]) ^ 1];
		}
	}
	b.Seq[count] = '\0';
	for (int i = 0; i != 256; ++i)
	{
		b.name[i] = name[i];
		if (name[i] == '\0') break;
	}
}

bool sequence::readNext(FILE * ifile)
{
	if (ifile == nullptr)
	{
		throw - 102;
		return false;
	}
	char buffer[BUFFER_SIZE];
	fscanf(ifile, "%s", buffer);
	if (feof(ifile)) return false;
	if (buffer[0] != '@')
	{
		throw - 112;//file content error
		return false;
	}
	int i;
	for (i = 0; i != 256; ++i)
	{
		name[i] = buffer[i + 1];
		if (buffer[i + 1] == '\0')
		{
			break;
		}
	}
	fscanf(ifile, "%s", buffer);
	for (i = 0; i != MAX_SEQUENCE_LENGTH; ++i)
	{
		if (buffer[i] == '\0') break;
		buffer[i] &= 0xDF;
		if (
			buffer[i] != 'A' &&
			buffer[i] != 'T' &&
			buffer[i] != 'G' &&
			buffer[i] != 'C' &&
			buffer[i] != 'N'
			)
		{
			throw - 112;
			return false;
		}
		Seq[i] = buffer[i];
	}
	Seq[i] = buffer[i];
	length = i;
	fscanf(ifile, "%s", buffer);
	fscanf(ifile, "%s", buffer);
	return true;
}

FILE* sequence::writeHead(const char *fname, index& theIndex)
{
	if (fname == nullptr)
	{
		return false;
	}
	FILE * file;
	file = fopen(fname, "w");
	if (file == nullptr)
	{
		throw - 102;
		return false;
	}
	fprintf(file, "@SQ\tSN:%s\tLN:%u", theIndex.name, theIndex.getLengthOfGenome());
	outputFile = file;
	return file;
}

inline int min(int a, int b)
{
	return a > b ? b : a;
}

bool sequence::write(FILE* ofile, index& theIndex)
{
	if (ofile == nullptr)
	{
		throw - 102;
		return false;
	}
	//char theQuilities[MAX_SEQUENCE_LENGTH + 1];//just some twos
	//for (int i = 0; i != length; ++i)
	//{
	//	theQuilities[i] = '2';
	//}
	//theQuilities[length] = '\0';
	All2s[length] = '\0';
	fprintf(ofile, "\n%s\t%d\t%s\t%u\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tPG:Z:SNAP\tRG:Z:FASTQ", name, isReversed ? 16 : 0, theIndex.name, Position, 60, CIGAR, "*", 0, 0, Seq, All2s);
	if (Position != 0)
	{
		fprintf(ofile, "\tNM:i:%u", editDis);
	}
	All2s[length] = '2';
	return true;
}

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

sequence reversed;
bool sequence::align(index & theIndex, sequence* &bestS)
{
	++SequenceAligned;
	size_t obest, osecond, opos, rbest, rsecond, rpos;
	int rreturn, oreturn;
	oreturn = singleAlign(theIndex, obest, osecond, opos);
	rc(reversed);
	rreturn = reversed.singleAlign(theIndex, rbest, rsecond, rpos);
	bestS = this;
	if (oreturn == -1 && rreturn == -1)
	{
		if (obest > rbest) bestS = &reversed;
		bestS->Position = 0;
	}
	else
	{
		if (oreturn == -1) bestS = &reversed;
		else if (rreturn != -1)
		{
			if (rbest < obest) bestS = &reversed;
		}
	}
	//if (doWrite)
	//{
	//	//while (WthreadLock);
	//	WthreadLock = true;
	//	if (outputFile == nullptr)
	//	{
	//		char *samFileName;
	//		if (samName[0] == '\0')
	//		{
	//			const char *samFileExtension = ".sam";
	//			size_t samFileNameLength = (size_t)(strlen(theIndex.name) + strlen(samFileExtension) + 1);
	//			samFileName = new char[samFileNameLength];
	//			_snprintf(samFileName, samFileNameLength, "%s%s", theIndex.name, samFileExtension);
	//		}
	//		else
	//		{
	//			const char *samFileExtension = ".sam";
	//			size_t samFileNameLength = (size_t)(strlen(samName) + strlen(samFileExtension) + 1);
	//			samFileName = new char[samFileNameLength];
	//			_snprintf(samFileName, samFileNameLength, "%s%s", samName, samFileExtension);
	//		}
	//		bestS->writeHead(samFileName, theIndex);
	//		bestS->write(outputFile, theIndex);
	//	}
	//	else
	//	{
	//		bestS->write(outputFile, theIndex);
	//	}
	//	WthreadLock = false;
	//}
	if (oreturn == -1 && rreturn == -1) return false;
	return true;
}

void sequence::writeResult(index& theIndex)
{
	//while (WthreadLock);
	if (outputFile == nullptr)
	{
		char *samFileName;
		if (samName[0] == '\0')
		{
			const char *samFileExtension = ".sam";
			size_t samFileNameLength = (size_t)(strlen(theIndex.name) + strlen(samFileExtension) + 1);
			samFileName = new char[samFileNameLength];
			_snprintf(samFileName, samFileNameLength, "%s%s", theIndex.name, samFileExtension);
		}
		else
		{
			const char *samFileExtension = ".sam";
			size_t samFileNameLength = (size_t)(strlen(samName) + strlen(samFileExtension) + 1);
			samFileName = new char[samFileNameLength];
			_snprintf(samFileName, samFileNameLength, "%s%s", samName, samFileExtension);
		}
		writeHead(samFileName, theIndex);
		write(outputFile, theIndex);
	}
	else
	{
		write(outputFile, theIndex);
	}
}

int sequence::singleAlign(index& theIndex, size_t &best, size_t &second, size_t &pos)
{
	size_t InitSeedsCount = 0, SeedsCount = 0, SeedsTried = 0;
	RSeed theSeedsToTry[MAX_SEED_NUM];
	size_t candidates[MAX_SEED_NUM*MAX_HIT_NUM], cn = 0, tried = 0, nHits;
	size_t * Hits;
	if (length < (size_t)seed::length)
	{
		throw - 114;//Read sequence is too short!
		return -1;
	}
	best = MAX_EDD;
	second = MAX_EDD;
	size_t round = 1;
	while (SeedsCount < SeedsToTry)
	{
		//gain seeds
		for (size_t i = round == 1 ? 0 : seed::length / round; i < length; i += seed::length)
		{
			if (SeedsCount >= SeedsToTry) break;
			try
			{
				theSeedsToTry[SeedsCount].acquireSeq(Seq + i);
				theSeedsToTry[SeedsCount].positionInSeq = i;
				++SeedsCount;
			}
			catch (int err)
			{
				if (err == -201)
				{
					break;
				}
				else
				{
					throw err;
					return -1;
				}
			}
		}
		if (round == 1)
		{
			InitSeedsCount = SeedsCount;
		}
		//gain candidates
		for (; SeedsTried < SeedsCount; ++SeedsTried)
		{
			nHits = theIndex.lookup(theSeedsToTry[SeedsTried].toNum(), Hits);
			//if (Hits == nullptr)
			//{
			//	free(Hits);
			//	throw - 115;//unknown error
			//	return -1;
			//}
			if (nHits == 0)
			{
				continue;
			}
			if (nHits > hMax)
			{
				continue;
			}
			size_t bpos;
			for (size_t j = 0; j != nHits; ++j)
			{
				//check if this position occured
				bool occured = false;
				bpos = (size_t)((theIndex.Index[Hits[j]].value - theSeedsToTry[SeedsTried].positionInSeq) / SAME_POS_BLOCK);
				for (size_t ci = 0; ci != cn; ++ci)
				{
					if (candidates[ci] / SAME_POS_BLOCK == bpos)
					{
						occured = true;
						break;
					}
				}
				if (occured)
				{
					continue;
				}
				candidates[cn++] = (size_t)(theIndex.Index[Hits[j]].value - theSeedsToTry[SeedsTried].positionInSeq);
			}
			free(Hits);
		}
		//calculate the edit distance
		for (; tried < cn; ++tried)
		{
			int edd = computeEditDistance(Seq, length, theIndex.getTheGenome() + candidates[tried] - 1, length, best);
			if (edd<0) continue;
			if (best>edd)
			{
				Position = candidates[tried];
				best = edd;
			}
			else if (second > edd)
			{
				second = edd;
			}
		}
		if (best <= dMax)
		{
			if (best + 2 <= InitSeedsCount)
			{
				if (best + 2 <= second)
				{
					editDis = best;
					return 0;
				}
				editDis = best;
				return 1;
			}
		}
	}
	if (best <= dMax)
	{
		if (best + 2 <= second)
		{
			editDis = best;
			return 0;
		}
		editDis = best;
		return 1;
	}
	editDis = best;
	return -1;
}

short L[MAX_K + 1][2 * MAX_K + 1];// TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
int sequence::computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k)
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
		uint64 x = *((uint64*)p) ^ *((uint64*)t);
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
					uint64 x = *((uint64*)p) ^ *((uint64*)t);
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

size_t sequence::SeedsToTry = 25;
size_t sequence::dMax = 8;
size_t sequence::c = 2;//confidence threshold
size_t sequence::hMax = 250;
char sequence::samName[256] = {};
bool sequence::doWrite = 1;
FILE* sequence::outputFile = nullptr;
size_t sequence::SequenceAligned = 0;