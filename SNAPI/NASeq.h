#pragma once

#include "base.h"
#include <stdio.h>

class NASeq
{
public:
	base* seq;
	unsigned LengthInByte;
	bool isSkipped;
	static char DNAcharset[16];
	static char RNAcharset[16];
public:
	bool isRNA;
	unsigned length;
	NASeq& reverse();
	NASeq& operator=(const NASeq&);
	NASeq& operator~();
	base operator[](unsigned) const;
	NASeq& append(const NASeq&);
	NASeq& append(const char *, unsigned = 0);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	char * toString(unsigned pos = 0, unsigned length = 0);//generate the at(u)gc string. while length=0, generate frome pos to the end of the sequence;
	//constructors
	NASeq();
	NASeq(const char*, unsigned=0);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	NASeq(const NASeq&);//equals to operator=
	NASeq(FILE*, bool isFQ=false);//receive fa file by default
	~NASeq();
};

