#include "seed.h"

unsigned seed::seedLen = 20;

bool seed::isASeed(const NASeq & S)
{
	if (S.length < seedLen) return false;
	for (unsigned i = 0; i != S.length; ++i)
	{
		if (NASeq::LegibleSeed[S.seq[i]] == 0) return false;
	}
	return true;
}

seed::seed(unsigned long long thevalue)
:value(thevalue)
{
}
//
//seed::seed(const NASeq& S)
//:value(0)
//{
//	for (unsigned i = 0; i < seedLen; i++)
//	{
//		if (i>S.length) break;
//		unsigned long long encodedBase = S[i].val>>4;
//		if (NASeq::DNAcharset[encodedBase] == 'V' || NASeq::DNAcharset[encodedBase] == 'N') throw - 1;
//		encodedBase &= 0x3;
//		value |= encodedBase << ((seedLen - i - 1) * 2);
//		//reversed |= (encodedBase ^ 0x3) << (i * 2);
//	}
//}

seed::~seed()
{
}

//seed& seed::reverse()
//{
//	forward += reversed;
//	reversed = forward - reversed;
//	forward -= reversed;
//	return *this;
//}
//seed seed :: operator~()
//{
//	return seed(reversed, forward);
//}
