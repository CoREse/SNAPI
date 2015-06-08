#pragma once

#include "base.h"
#include <stdio.h>
#include <map>

//NASeq copies(only copies) share the same dynamic allocated memory, that's to say, you change one of them, you change all of them. That's awkward, it was changed to this way to improve the efficiency
class NASeq
{
	static std::map <char*, unsigned> Copies;//key is the originalSeq, and value is the number of copies
	//unsigned long copies;//the number of this NASeq's copies, when it's 0, it's time to free the space
	char* OriginalSeq;//pointed to the original seq(for subseqs)

	bool isSubString;
	bool isQuickCopyable;

	static char buffer[10240];
public:
	char* seq;
	unsigned length;
	static int LegibleSeed[256];
	static char Opposite[256];
	static unsigned NANumber[256];

	static int computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k, short (*L)[63]);

	NASeq& reverse();
	NASeq& operator=(const NASeq&);//because we share the same dynamic space, I must change your copies attribute
	NASeq operator~();
	NASeq getSubSequence(int start_loc, unsigned length);// length must be above 0
	char operator[](unsigned) const;

	bool readNext(FILE*, char* =nullptr);

	//future work
	bool isASeed(unsigned start=0, unsigned SeedLength=0) const;//length=0 means the whole length
	unsigned long long getSeed(unsigned start=0, unsigned SeedLength=0) const;//length=0 means the whole length, shall throw -201 when reach the end

	//append is now not safe for copying
	//NASeq& append(const NASeq&);
	//NASeq& append(const char *, unsigned = 0);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	//bool toString(char * str, unsigned pos = 0, unsigned length = 0);//generate the at(u)gc string. while length=0, generate frome pos to the end of the sequence;
	//constructors
	NASeq();
	NASeq(const char*, unsigned=0, bool isQuickCopyable=false);//if the string contains letters other than atugcnATUGCN then it will cause unexpected errors, we won't check it, but the letters in the string shall be all atugcATUGCnN
	NASeq(const NASeq&, bool isQuickCopyable=false);//equals to operator=
	NASeq(FILE*, bool isFQ=false, char* seqName=nullptr);//receive fa file by default, non-thread save
	~NASeq();

	bool saveToFile(FILE*);
	bool loadFromFile(FILE*);
};