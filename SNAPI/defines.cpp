#include "defines.h"
#include <stdio.h>

//the errors strings, the index=-error-101.
char errors[15][128] =
{
	"Some letter is not a Base but is treated as a Base!",
	"File pointer is null!",
	"The genome reference's length is too short!",
	"The name read from fa file is too long!",
	"Failed to insert into Index!",
	"Failed to malloc!",
	"Something wrong with the data in the Index!"
	"Abnormal islast flag!",
	"An error occured while allocating space!",
	"Failed to write Index into file!",
	"Failed to read Index from file!",
	"File content error!",
	"The K in compute edit distance is too big!",
	"The read sequence is too short!",
	"Unknown error!"
};
//May cause error if err is not a error
void printError(int err)
{
	printf("%s",errors[-err - 101]);
}