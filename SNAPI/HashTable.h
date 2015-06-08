#include <vector>
//HashTable using square detection
class HashTable
{
public:
	static float OverflowPara, MainPara;//nOverflowTable/NumberOfEntries and nMainTable/NumberOfEntries
	static unsigned long long UnusedKey;
	static unsigned UnusedLocation;
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
		unsigned location1;
		unsigned location2;
	};
	struct OverflowEntry
	{
		unsigned long long key;
		std::vector<unsigned> locations;
	};
	class Result
	{
		Entry* Main;
		OverflowEntry* Overflow;
	public:
		Result(Entry* MainResult = nullptr, OverflowEntry* OverflowResult = nullptr);
		unsigned size() const
		{
			if (Main == nullptr) return 0;
			if (Main->location1 == UnusedLocation) return 0;
			if (Main->location2 == UnusedLocation) return 1;
			if (Overflow == nullptr) return 2;
			return Overflow->locations.size() + 2;
		}
		unsigned& operator[](unsigned index);
	};
	~HashTable();
	HashTable(unsigned NumberOfEntries = 0);
	inline unsigned getnMainTable() const { return nMainTable; }
	inline unsigned getnUsedMainTable() const { return nUsedMainTable; }
	inline unsigned getnOverflowTable() const { return nOverflowTable; }
	inline unsigned getnUsedOverflowTable() const { return nUsedOverflowTable; }

	bool insert(unsigned long long key, unsigned location);
	Entry* lookupMainTable(unsigned long long key);
	OverflowEntry * lookupOverflowTable(unsigned long long key);
	//OverflowEntry lookupIntel(unsigned long long) const;
	Result lookup(unsigned long long key);

	bool saveToFile(FILE*);
	bool loadFromFile(FILE*);
private:
	Entry *MainTable;
	OverflowEntry * OverflowTable;
	unsigned nMainTable, nOverflowTable, nUsedMainTable, nUsedOverflowTable;
	friend class Indexer;
};
