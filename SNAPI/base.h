#pragma once

//struct __base
//{
//	unsigned char low : 4;
//	unsigned char high : 4;
//};
//
//union _base
//{
//	__base base;
//	unsigned char val;
//};

struct base
{
	unsigned char val;
	friend class sequence;
	base();
	base(char);
	base(const char*);
	base(int);
	base & operator =(const char *);
	base & operator =(const base &);
	base& operator ~();
	base & reverse();
	~base();
};

