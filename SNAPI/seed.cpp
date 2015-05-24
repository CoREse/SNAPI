#include "seed.h"

unsigned seed::seedLen = 20;

bool seed::isASeed(const NASeq & S)
{
	if (S.length < seedLen) return false;
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
