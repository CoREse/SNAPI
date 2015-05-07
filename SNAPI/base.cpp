#include "base.h"


base::base()
:val(0)
{
}

base::base(char _val)
{
	_val &= 0xDF;//upper case
	switch (_val)
	{
	case 'A':
		val = 0x80;
		break;
	case 'T':
	case 'U':
		val = 0xB0;
		break;
	case 'G':
		val = 0x90;
		break;
	case 'C':
		val = 0xa0;
		break;
	case 'N':
		val = 0xc0;
		break;
	default:
		throw 1;
		break;
	}
}

base::base(const char * seq)
{
	switch ((*seq) & 0xDF)
	{
	case 'A':
		val = 0x80;
		break;
	case 'T':
	case 'U':
		val = 0xB0;
		break;
	case 'G':
		val = 0x90;
		break;
	case 'C':
		val = 0xa0;
		break;
	case 'N':
		val = 0xc0;
		break;
	default:
		throw 1;
		break;
	}
	switch ((*(seq + 1)) & 0xDF)
	{
	case 'A':
		val |= 0x08;
		break;
	case 'T':
	case 'U':
		val |= 0x0B;
		break;
	case 'G':
		val |= 0x09;
		break;
	case 'C':
		val |= 0x0a;
		break;
	case 'N':
		val = 0x0c;
		break;
	default:
		return;
		break;
	}
}

base::~base()
{
}

base& base::reverse()
{
	val ^= 0x33;
	return *this;
}

base & base::operator =(const char * sval)
{
	switch ((*sval) & 0xDF)
	{
	case 'A':
		val = 0x80;
		break;
	case 'T':
	case 'U':
		val = 0xB0;
		break;
	case 'G':
		val = 0x90;
		break;
	case 'C':
		val = 0xa0;
		break;
	case 'N':
		val = 0xc0;
		break;
	default:
		throw 1;
		break;
	}
	switch ((*(sval + 1)) & 0xDF)
	{
	case 'A':
		val |= 0x08;
		break;
	case 'T':
	case 'U':
		val |= 0x0B;
		break;
	case 'G':
		val |= 0x09;
		break;
	case 'C':
		val |= 0x0a;
		break;
	case 'N':
		val |= 0x0c;
		break;
	default:
		return *this;
		break;
	}
	return *this;
}

base & base::operator =(const base & abase)
{
	val = abase.val;
	return *this;
}

base::base(int aint)
:val(aint)
{
}

base& base::operator~()
{
	return reverse();
}