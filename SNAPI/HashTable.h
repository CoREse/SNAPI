#include <vector>
class HashTable
{
Entry *MainTable;
OverfolwEntry * OverflowTable;
unsigned long long nMainTable, nOverflowTable, nUsedMainTable, nUsedOverflowTable;
public:
static double OverflowPara;
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
HashTable(unsigned long long NumberOfEntries=0);
inline unsigned long long getnMainTable() const (return nMainTable);
inline unsigned long long getnUsedMainTable() const (return nUsedMainTable);
inline unsigned long long getnOverflowTable() const (return nOverflowTable);
inline unsigned long long getnUsedOverflowTable() const (return nUsedOverflowTable);

bool insert(unsigned long long key, unsigned long long location);
Entry* lookupMainTable(unsigned long long key);
OverflowEntry * lookupOverflowTable(unsigned long long key);
OverflowEntry lookupIntel(unsigned long long) const;
void lookup(unsigned long long key, Entry* MainResult, Overflow* OverflowResult);
};
