#include <vector>
//HashTable using square detection
class HashTable
{
public:
	static double OverflowPara;//nOverflowTable/nMainTable
	static unsigned long long UnusedKey;
	static unsigned long long Hash(unsigned long long key)
	{
		key ^= key >> 33;
		key *= 0xff51afd7ed558ccd;
		key ^= key >> 33;
		key *= 0xc4ceb9fe1a85ec53;
		key ^= key >> 33;
		return key;
	}
	struct Entry{
		unsigned long long key;
		unsigned long long location;
	};
	struct OverflowEntry
	{
		unsigned long long key;
		std::vector<unsigned long long> locations;
	};
	~HashTable();
	HashTable(unsigned long long NumberOfEntries = 0);
	inline unsigned long long getnMainTable() const { return nMainTable; }
	inline unsigned long long getnUsedMainTable() const { return nUsedMainTable; }
	inline unsigned long long getnOverflowTable() const { return nOverflowTable; }
	inline unsigned long long getnUsedOverflowTable() const { return nUsedOverflowTable; }

	bool insert(unsigned long long key, unsigned long long location);
	Entry* lookupMainTable(unsigned long long key);
	OverflowEntry * lookupOverflowTable(unsigned long long key);
	//OverflowEntry lookupIntel(unsigned long long) const;
	void lookup(unsigned long long key, Entry* MainResult, OverflowEntry* OverflowResult);
private:
	Entry *MainTable;
	OverflowEntry * OverflowTable;
	unsigned long long nMainTable, nOverflowTable, nUsedMainTable, nUsedOverflowTable;
};
