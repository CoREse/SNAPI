#include "NASeq.h"
#include <stdlib.h>
#include <string.h>
#include <cassert>

#define MAX_SEQUENCE_LENGTH 3000000000

char buffer[10240];
//char DNAcharset[17] = "VVVVVVVVAGCTNNNN";
char NASeq::DNAcharset[16] =
{
	'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V',//V means vacancy
	'A', 'G', 'C', 'T', 'N', 'N', 'N', 'N'
};
char NASeq::RNAcharset[16] =
{
	'V', 'V', 'V', 'V', 'V', 'V', 'V', 'V',//V means vacancy
	'A', 'G', 'C', 'U', 'N', 'N', 'N', 'N'
};

NASeq::NASeq()
:isRNA(false), isSkipped(false), LengthInByte(0), length(0), seq(nullptr), copies(1), originalSeq(nullptr)
{
}

NASeq::NASeq(const char * sseq, unsigned slength)
: isRNA(false), isSkipped(false), copies(1), originalSeq(nullptr)
{
	length = slength == 0 ? strlen(sseq) : slength;
	if (length == 0)
	{
		LengthInByte = 0;
		seq = nullptr;
		return;
	}
	LengthInByte = length / 2 + length % 2;
	seq = (base*)malloc(LengthInByte*sizeof(base));
	for (unsigned i = 0; i < LengthInByte; ++i)
	{
		seq[i] = sseq + i * 2;
	}
	originalSeq = seq;
}

NASeq::NASeq(NASeq& b)
{
	*this = b;
}

NASeq::NASeq(FILE * ifile, bool isFQ)
:copies(1), originalSeq(nullptr)
{
	if (isFQ)
	{
		if (ifile == nullptr)
		{
			throw - 102;
			return;
		}
		fscanf(ifile, "%s", buffer);
		if (feof(ifile)) return;
		if (buffer[0] != '@')
		{
			throw - 112;//file content error
			return;
		}
		int i;
		fscanf(ifile, "%s", buffer);
		NASeq(buffer);
		fscanf(ifile, "%s", buffer);
		fscanf(ifile, "%s", buffer);
	}
	else
	{
		if (ifile == nullptr)
		{
			throw - 102;//file pointer is null
			return;
		}
		//allocate space
		length = 0;
		LengthInByte = 0;
		fseek(ifile, 0, SEEK_END);
		seq = new base[ftell(ifile) / 2 + 1];//the sizeof(base) is 1. The fa file always provide a very long sequence, in this case, we don't mind to allocate a little more space for restoring it
		if (seq == nullptr)
		{
			fclose(ifile);
			throw - 106;//malloc·ÖÅä¿Õ¼äÊ§°Ü
			return;
		}
		rewind(ifile);
		//read the sequence
		fscanf(ifile, "%s", buffer);//the first line is name of the sequence, which in this class we don't care
		while (true)
		{
			fscanf(ifile, "%s", buffer);
			if (feof(ifile)) return;
			//append(buffer);//very slow
			unsigned slen = strlen(buffer);
			if (slen == 0)
				continue;
			if (isSkipped == (length % 2 == 0))
				//Ä©Î²ÓÐ¿Õ
			{
				seq[LengthInByte - 1].val &= 0xf0;
				switch ((*buffer) & 0xDF)
				{
				case 'A':
					seq[LengthInByte - 1].val += 0x08;
					break;
				case 'T':
				case 'U':
					seq[LengthInByte - 1].val += 0x0B;
					break;
				case 'G':
					seq[LengthInByte - 1].val += 0x09;
					break;
				case 'C':
					seq[LengthInByte - 1].val += 0x0a;
					break;
				case 'N':
					seq[LengthInByte - 1].val += 0x0c;
					break;
				default:
					throw 1;
					break;
				}
				for (unsigned i = 0; i != (slen - 1) / 2; ++i)
				{
					seq[LengthInByte + i] = buffer + 1 + i * 2;
				}
				if (slen % 2 == 0)
				{
					seq[LengthInByte + slen / 2 - 1] = buffer[slen - 1];
				}
				LengthInByte += slen / 2;
			}
			else
				//Ä©Î²ÎÞ¿Õ
			{
				for (unsigned i = 0; i != slen / 2 + slen % 2; ++i)
				{
					seq[LengthInByte + i] = buffer + 2 * i;
				}
				LengthInByte += slen / 2 + slen % 2;
			}
			length += slen;
		}
	}
	originalSeq = seq;
}

NASeq::~NASeq()
{
	if (--copies == 0)	delete originalSeq;
}

NASeq& NASeq::operator=(NASeq & b)
{
	++b.copies;
	originalSeq = b.originalSeq;
	isRNA = b.isRNA;
	isSkipped = b.isSkipped;
	length = b.length;
	LengthInByte = b.LengthInByte;
	seq = b.seq;
	copies = b.copies;
	return *this;
}

/*the old version
NASeq& NASeq::operator=(NASeq & b)
{
	++b.copies;
	isRNA = b.isRNA;
	isSkipped = b.isSkipped;
	length = b.length;
	LengthInByte = b.LengthInByte;
	seq = (base*)malloc(LengthInByte*sizeof(base));
	for (unsigned i = 0; i != LengthInByte; ++i)
	{
		seq[i] = b.seq[i];
	}
	copies = b.copies;
	return *this;
}
*/

NASeq NASeq::operator~()
{
	return reverse();
}

NASeq& NASeq::reverse()
{
	for (unsigned i = 0; i != LengthInByte; ++i)
	{
		~seq[i];
	}
	return *this;
}

/*
NASeq& NASeq::append(const NASeq& b)
{
	if (b.length == 0)
		return *this;
	if (isSkipped == (length % 2 == 0))
		//Ä©Î²ÓÐ¿Õ
	{
		seq = (base*)realloc(seq, LengthInByte + b.length / 2);//sizeof(base) is 1
		if (b.isSkipped)
		{
			seq[LengthInByte - 1].val &= 0xf0;
			seq[LengthInByte - 1].val |= b.seq[0].val & 0x0f;
			for (unsigned i = 1; i < b.LengthInByte; ++i)
			{
				seq[LengthInByte + i - 1] = b.seq[i];
			}
		}
		else
		{
			seq[LengthInByte - 1].val &= 0xf0;
			seq[LengthInByte - 1].val |= b.seq[0].val >> 4;
			for (unsigned i = 0; i != b.LengthInByte - 1; ++i)
			{
				seq[LengthInByte + i].val = b.seq[i].val << 4 + b.seq[i + 1].val >> 4;
			}
			if (b.length % 2 == 0)
			{
				seq[LengthInByte + b.LengthInByte - 1].val = b.seq[b.LengthInByte - 1].val << 4;
			}
		}
		LengthInByte += b.length / 2;
	}
	else
		//Ä©Î²ÎÞ¿Õ
	{
		seq = (base*)realloc(seq, LengthInByte + b.length / 2 + b.length % 2);//sizeof(base) is 1
		if (b.isSkipped)
		{
			for (unsigned i = 0; i != b.LengthInByte - 1; ++i)
			{
				seq[LengthInByte + i].val = b.seq[i].val << 4 | (b.seq[i + 1].val >> 4);
			}
			if (b.length % 2 == 1)
			{
				seq[LengthInByte + b.LengthInByte - 1].val = b.seq[LengthInByte - 1].val << 4;
			}
		}
		else
		{
			for (unsigned i = 0; i != b.LengthInByte; ++i)
			{
				seq[LengthInByte + i] = b.seq[i];
			}
		}
		LengthInByte += b.length / 2 + b.length % 2;
	}
	length += b.length;
	return *this;
}

NASeq& NASeq::append(const char * sseq, unsigned slen)
{
	if (slen == 0)
	{
		slen = strlen(sseq);
	}
	if (slen == 0)
		return *this;
	if (isSkipped == (length % 2 == 0))
		//Ä©Î²ÓÐ¿Õ
	{
		seq = (base*)realloc(seq, LengthInByte + slen / 2);//sizeof(base) is 1
		seq[LengthInByte - 1].val &= 0xf0;
		switch ((*sseq) & 0xDF)
		{
		case 'A':
			seq[LengthInByte - 1].val += 0x08;
			break;
		case 'T':
		case 'U':
			seq[LengthInByte - 1].val += 0x0B;
			break;
		case 'G':
			seq[LengthInByte - 1].val += 0x09;
			break;
		case 'C':
			seq[LengthInByte - 1].val += 0x0a;
			break;
		case 'N':
			seq[LengthInByte - 1].val += 0x0c;
			break;
		default:
			throw 1;
			break;
		}
		for (unsigned i = 0; i != (slen - 1) / 2; ++i)
		{
			seq[LengthInByte + i] = sseq + 1 + i * 2;
		}
		if (slen % 2 == 0)
		{
			seq[LengthInByte + slen / 2 - 1] = sseq[slen - 1];
		}
		LengthInByte += slen / 2;
	}
	else
		//Ä©Î²ÎÞ¿Õ
	{
		seq = (base*)realloc(seq, LengthInByte + slen / 2 + slen % 2);//sizeof(base) is 1
		for (unsigned i = 0; i != slen / 2 + slen % 2; ++i)
		{
			seq[LengthInByte + i] = sseq + 2 * i;
		}
		LengthInByte += slen / 2 + slen % 2;
	}
	length += slen;
	return *this;
}
*/

/*
char * NASeq::toString(unsigned pos, unsigned length)
{
	if (length == 0)
	{
		length = this->length - pos;
	}
	//#################################################
	//note: need futher job to handle the memory
	//####################################################
	char * str = new char[length + 1];
	if (isRNA)
	{
		if (isSkipped)
		{
			for (unsigned i = 0; i != length; ++i)
			{
				if (pos + i >= this->length) break;
				str[i] = RNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val & 0x0F : seq[(pos + i) / 2 + 1].val >> 4];
			}
		}
		else
		{
			for (unsigned i = 0; i != length; ++i)
			{
				if (pos + i >= this->length) break;
				str[i] = RNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val >> 4 : seq[(pos + i) / 2].val & 0x0F];
			}
		}
	}
	else
	{
		if (isSkipped)
		{
			for (unsigned i = 0; i != length; ++i)
			{
				if (pos + i >= this->length) break;
				str[i] = DNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val & 0x0F : seq[(pos + i) / 2 + 1].val >> 4];
			}
		}
		else
		{
			for (unsigned i = 0; i != length; ++i)
			{
				if (pos + i >= this->length) break;
				str[i] = DNAcharset[(pos + i) % 2 == 0 ? seq[(pos + i) / 2].val >> 4 : seq[(pos + i) / 2].val & 0x0F];
			}
		}
	}
	str[length] = '\0';
	return str;
}
*/

base NASeq::operator[] (unsigned i) const
{
	if (isSkipped)
	{
		return base(seq[i / 2 + i % 2].val << i % 2 == 0 ? 4 : 0);
	}
	else
	{
		return base(seq[i / 2].val << i % 2 == 0 ? 0 : 4);
	}
}


NASeq NASeq::getSubSequence(unsigned start_loc, unsigned length)
{
	assert(length != 0);
	NASeq sub;
	sub.copies = ++copies;
	sub.isRNA = isRNA;
	sub.length = length;
	sub.originalSeq = originalSeq;
	if (isSkipped)
	{
		sub.seq = seq + length / 2 + length % 2;
		sub.isSkipped = start_loc % 2 == 0;
		sub.LengthInByte = sub.isSkipped ? 1 : 0 + length / 2 + sub.isSkipped ? 0 : length % 2;
	}
	else
	{
		sub.seq = seq + start_loc / 2;
		sub.isSkipped = start_loc % 2;
		sub.LengthInByte = sub.isSkipped ? 1 : 0 + length / 2 + sub.isSkipped ? 0 : length % 2;
	}
	sub.copies = copies;
	return *this;
}


bool NASeq::saveToFile(FILE* file)
{
	fprintf(file, "%d%d%u%u", isSkipped, isRNA, length, LengthInByte);
	if (fwrite(seq, sizeof(base), LengthInByte, file) != LengthInByte)
	{
		fclose(file);
		throw - 110;//Ð´ÈëNASeqÊ§°Ü£¡
		return false;
	}
	return true;
}

bool NASeq::loadFromFile(FILE* file)
{
	if (seq != nullptr) return false;
	fscanf(file, "%d%d%u%u", &isSkipped, &isRNA, &length, &LengthInByte);
	copies = 0;
	seq = (base*)malloc(LengthInByte);
	originalSeq = seq;
	if (fread(seq, sizeof(base), LengthInByte, file) != LengthInByte)
	{
		fclose(file);
		throw - 110;//Ð´ÈëNASeqÊ§°Ü£¡
		return false;
	}
	return true;
}