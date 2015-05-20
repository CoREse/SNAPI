#include "HashTable.h"
#include <stdlib.h>

HashTable::HashTable(unsigned long long NumberOfEntries)
:nMainTable(NumberOfEntries), nOverflowTable(NumberOfEntries*OverflowPara), nUsedMainTable(0), nUsedOverflowTable(0)
{
if (NumberOfEntries!=0)
{
MainTable=(Entry *)malloc(NumberOfEntries*sizeof(Entry));
OverflowTable=(OverflowEntry *) malloc(nOverflowTable*sizeof(OverflowEntry));
}
}
