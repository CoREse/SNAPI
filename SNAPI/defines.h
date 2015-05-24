/*
defines of the errors:
-101: some letter is not a Base but is treated as a Base
-102: file pointer is null
-103: the genome reference's length is too short
-104: the name from fa file is too long
-105: fail to insert into Index
-106: failed to malloc
-107: something wrong with the data in the Index
-108: abnormal islast flag
-109: an error occured while allocating space
-110: failed to write Index into file
-111: failed to read Index from file
-112: file content error
-113: the K in compute edit distance is too big
-114: read sequence is too short
-115: unknown error

*/

void printError(int);//receive error number and print the error.
//
//
#define MAX_SEED_LENGTH 32
#define uint64 unsigned long long
#define BUFFER_SIZE 1000
#define PROBE_DEPTH 4
#define MALLOC_BLOCK 1024
#define MAX_SEQUENCE_LENGTH 1024
//if pos/this is the same, then see them as a same position
#define SAME_POS_BLOCK 32
#define MAX_K 31
#define MAX_SEED_NUM 128
#define MAX_HIT_NUM 2048
#define MAX_EDD 1024
#define MAX_SAVE_SIZE 50000000u