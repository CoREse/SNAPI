#pragma once
#include "defines.h"
#include "NASeq.h"

struct seed
{
	static int seedLen;
	static bool isASeed(const NASeq &);
	//inline seed(const char * =nullptr);
	inline seed(const NASeq&);
	inline seed(unsigned long long forward, unsigned long long reversed);
	inline ~seed();
	inline seed& reverse();
	inline seed operator~();

	inline unsigned long long getBases() const
	{
		return forward;
	}

	inline unsigned long long getReversed() const
	{
		return reversed;
	}

	inline unsigned getLowBits() const
	{
		return forward;
	}
	inline unsigned getHighBits() const
	{
		return forward >> 32;
	}
	
private:
	unsigned long long forward;
	unsigned long long reversed;
};
