#include "seed.h"

int seed::seedLen = 20;

bool seed::isASeed(const NASeq & S)
{
	for (unsigned i = 0; i != S.length; ++i)
	{
		if (S[i].val >> 6 != 2) return false;
	}
	return true;
}

seed::seed(const NASeq& S)
:forward(0), reversed(0)
{
	for (unsigned i = 0; i < seedLen; i++)
	{
		unsigned long long encodedBase = S[i].val>>4;
		if (NASeq::DNAcharset[encodedBase] == 'V' || NASeq::DNAcharset[encodedBase] == 'N') throw - 1;
		encodedBase &= 0x3;
		forward |= encodedBase << ((seedLen - i - 1) * 2);
		reversed |= (encodedBase ^ 0x3) << (i * 2);
	}
}

seed::seed(unsigned long long f,unsigned long long r)
:forward(f), reversed(r)
{
}

seed::~seed()
{
}

seed& seed::reverse()
{
	forward += reversed;
	reversed = forward - reversed;
	forward -= reversed;
	return *this;
}
seed seed :: operator~()
{
	return seed(reversed, forward);
}

char bases[4] = { 'A', 'T', 'G', 'C' };

unsigned long long seed::toNum()//A=0,T=1,G=2,C=4,四进制数
{
	long long ret = 0;
	for (int i = 0; i <= seed::length; ++i)
	{
		if (Seq[i] != 'A'&&Seq[i] != 'T'&&Seq[i] != 'G'&&Seq[i] != 'C')
		{
			if (Seq[i] == 'N')
			{
				ret >>= 2;
				continue;
			}
			return ret;
		}
		ret += seed::convertBase2Num(Seq[i])*four[seed::length - i];
	}
	return ret;
}

unsigned long long const seed::four[32] =
{
	1,
	4,
	16,
	64,
	256,
	1024,
	4096,
	4096 * 4,
	4096 * 16,
	4096 * 64,
	4096 * 256,
	4096 * 1024,
	4096 * 4096,
	4096 * 4096 * 4,
	4096 * 4096 * 16,
	4096 * 4096 * 64,
	4096ui64 * 4096 * 256,
	4096ui64 * 4096 * 1024,
	4096ui64 * 4096 * 4096,
	4096ui64 * 4096 * 4096 * 4,
	4096ui64 * 4096 * 4096 * 16,
	4096ui64 * 4096 * 4096 * 64,
	4096ui64 * 4096 * 4096 * 256,
	4096ui64 * 4096 * 4096 * 1024,
	4096ui64 * 4096 * 4096 * 4096,
	4096ui64 * 4096 * 4096 * 4096 * 4,
	4096ui64 * 4096 * 4096 * 4096 * 16,
	4096ui64 * 4096 * 4096 * 4096 * 64,
	4096ui64 * 4096 * 4096 * 4096 * 256,
	4096ui64 * 4096 * 4096 * 4096 * 1024,
	4096ui64 * 4096 * 4096 * 4096 * 4096,
	4096ui64 * 4096 * 4096 * 4096 * 4096 * 4
};

void seed::acquireSeq(const char* S)
{
	if (S != nullptr)
	{
		int numberofbreak = 0;
		for (int i = 0; i < seed::length; ++i)
		{
			if (S[i+numberofbreak] == '\n')
			{
				++numberofbreak;
				--i;
				continue;
			}
			if (S[i + numberofbreak] == '\0') throw - 201;//reached the end of input
			Seq[i] = S[i + numberofbreak] & 0xDF;//只保留大写的基因字母
		}
	}
}

int seed::length = 20;//initialization for length