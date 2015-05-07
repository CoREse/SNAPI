#pragma once
#include "seed.h"
class RSeed :
	public seed
{
public:
	RSeed();
	~RSeed();
	size_t positionInSeq;//starts at 0
};

