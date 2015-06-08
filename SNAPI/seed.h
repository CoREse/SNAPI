#pragma once
#include "defines.h"
#include "NASeq.h"

struct seed
{
	static unsigned seedLen;
	static bool isASeed(const NASeq &);
	//inline seed(const char * =nullptr);
	//seed(const NASeq&);
	seed(const unsigned long long thevalue = 0);
	//seed(unsigned long long forward = 0, unsigned long long reversed = 0);
	~seed();
	//inline seed& reverse();
	//inline seed operator~();

	inline unsigned long long getBases() const
	{
		return value;
	}

	//inline unsigned long long getReversed() const
	//{
	//	return reversed;
	//}

	inline unsigned getLowBits() const
	{
		return value;
	}
	inline unsigned getHighBits() const
	{
		return value >> 32;
	}

	operator unsigned long long() const
	{
		return value;
	}
	
private:
	unsigned long long value;
//	unsigned long long forward;
//	unsigned long long reversed;
};
