#pragma once

#include "base.h"
#include <stdio.h>
#include <map>

//NASeq copies(only copies) share the same dynamic allocated memory, that's to say, you change one of them, you change all of them. That's awkward, it was changed to this way to improve the efficiency
class NASeq
{
	static std::map <base*, unsigned> copies;//key is the originalSeq, and value is the number of copies
	//unsigned long copies;//the number of this NASeq's copies, when it's 0, it's time to free the space
	base* originalSeq;//pointed to the original seq(for subseqs)
	bool isSubString;

	unsigned long long get64bit(unsigned index)const;
	unsigned long long get32bit(unsigned index)const;
public:
	base* seq;
	unsigned LengthInByte;
	bool isSkipped;
	static char DNAcharset[16];
	static char RNAcharset[16];

	static int computeEditDistance(const NASeq&, const NASeq&, int);

	bool isRNA;
	unsigned length;
	NASeq& reverse();
	NASeq& operator=(const NASeq&);//because we share the same dynamic space, I must change your copies attribute
	NASeq operator~();
	NASeq getSubSequence(int start_loc, unsigned length);// length must be above 0
	base operator[](unsigned) const;

	//future work
	//bool isASeed() const;
	//unsigned long long toSeed() const;

	//append is now not safe for copying
	//NASeq& append(const NASeq&);
	//NASeq& append(const char *, unsigned = 0);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	//char * toString(unsigned pos = 0, unsigned length = 0);//generate the at(u)gc string. while length=0, generate frome pos to the end of the sequence;
	//constructors
	NASeq();
	NASeq(const char*, unsigned=0);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	NASeq(const NASeq&);//equals to operator=
	NASeq(FILE*, bool isFQ=false);//receive fa file by default
	~NASeq();

	bool saveToFile(FILE*);
	bool loadFromFile(FILE*);
};