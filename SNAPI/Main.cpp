//#include <iostream>
#include "seed.h"
#include <stdio.h>
#include <time.h>
#include <vector>
#include "base.h"
#include "NASeq.h"
#include "Genome.h"
#include "Indexer.h"
#include "Aligner.h"
#include <omp.h>

using namespace std;

const char *SNAPI_VERSION = "0.01";
time_t start, stop;

unsigned threadNum=1;

void initiate()
{
	NASeq::LegibleSeed['A'] = 1;
	NASeq::LegibleSeed['T'] = 1;
	NASeq::LegibleSeed['G'] = 1;
	NASeq::LegibleSeed['C'] = 1;
	NASeq::LegibleSeed['a'] = 1;
	NASeq::LegibleSeed['t'] = 1;
	NASeq::LegibleSeed['g'] = 1;
	NASeq::LegibleSeed['c'] = 1;
	NASeq::NANumber['A'] = 0;
	NASeq::NANumber['T'] = 3;
	NASeq::NANumber['G'] = 1;
	NASeq::NANumber['C'] = 2;
	NASeq::NANumber['a'] = 0;
	NASeq::NANumber['t'] = 3;
	NASeq::NANumber['g'] = 1;
	NASeq::NANumber['c'] = 2;
	NASeq::Opposite['A'] = 'T';
	NASeq::Opposite['T'] = 'A';
	NASeq::Opposite['G'] = 'C';
	NASeq::Opposite['C'] = 'G';
	NASeq::Opposite['a'] = 't';
	NASeq::Opposite['t'] = 'a';
	NASeq::Opposite['g'] = 'c';
	NASeq::Opposite['c'] = 'g';
	Aligner::SeedsStartPositions[1] = 0.0;
	Aligner::SeedsStartPositions[2] = 0.25;
	Aligner::SeedsStartPositions[3] = 0.75;
	Aligner::SeedsStartPositions[4] = 0.125;
	Aligner::SeedsStartPositions[5] = 3.0 / 8.0;
	Aligner::SeedsStartPositions[6] = 5.0 / 8.0;
	Aligner::SeedsStartPositions[7] = 7.0 / 8.0;
	Aligner::SeedsStartPositions[8] = 1.0 / 16.0;
	Aligner::SeedsStartPositions[9] = 3.0 / 16.0;
	Aligner::SeedsStartPositions[10] = 5.0 / 16.0;
	Aligner::SeedsStartPositions[11] = 7.0 / 16.0;
	Aligner::SeedsStartPositions[12] = 9.0 / 16.0;
	Aligner::SeedsStartPositions[13] = 11.0 / 16.0;
	Aligner::SeedsStartPositions[14] = 13.0 / 16.0;
	Aligner::SeedsStartPositions[15] = 15.0 / 16.0;

}

static void usage()
{
	fprintf(stderr,
		"Usage: snap <command> [<options>]\n"
		"Commands:\n"
		"   index    build a genome index\n"
		"   align    align single-end reads\n"
		"Type a command without arguments to see its help.\n");
	exit(1);
}

void buildIndex(int argc, const char **argv);
void runAlign(int, const char **);

int main(int argc, const char **argv)
{
	initiate();
	//Indexer theIndexer;
	////theIndexer.readReference("G:\\");
	////theIndexer.buildIndex();

	//theIndexer.loadFromFile("F:\\save\\OriginalChr1\\save");
	////theIndexer.loadFromFile("F:\\save\\3rdCutEdition\\hMax100chr1\\cutedsave");
	//
	////theIndexer.throughlyCut(0, "F:\\save\\ThroughlyCut\\chr21\\cutedsave");

	//theIndexer.cutIndex(Indexer::hMax, "F:\\save\\3rdCutEdition\\hMax100chr1\\cutedsave");
	////theIndexer.saveToFile("F:\\save\\OriginalChr1\\save");
	////theIndexer.cutIndex(1000);
	////char astr[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	////uint64 bstr = *(uint64*)astr;

	//Aligner theAligner;
	//theAligner.singleAlign("E:\\Programs\\Bioinformatics\\data\\hg19\\chr1.fq","G:\\result.sam");
	////theAligner.singleAlign("G:\\chr21.fq", "G:\\result.sam");
	//return 0;
//	//unsigned index;
//	//short L[MAX_K + 1][2 * MAX_K + 1];
//	//NASeq aSeq("GGGGACTGTTGTGGGGTGGGGGGAGGGGGGAGGGATAGCATTGGGAGATATACCTAATGCTAGATGATGAGTTAGTGGGTGCAGCGCACCAGCATGGCAC");
//	//aSeq.reverse();
//	//unsigned UnscoredMostHitedLocation = 9656354;
//	//unsigned chrNum = Indexer::Reference->getChrNumber(UnscoredMostHitedLocation);
//	//cout<<NASeq::computeEditDistance(Indexer::Reference->chrs[chrNum].sequence.seq + UnscoredMostHitedLocation - Indexer::Reference->chrs[chrNum].start_location, aSeq.length, aSeq.seq, aSeq.length, 100, L);
//	//UnscoredMostHitedLocation = 23032888;
//	//chrNum = Indexer::Reference->getChrNumber(UnscoredMostHitedLocation);
//	//cout << endl << NASeq::computeEditDistance(Indexer::Reference->chrs[chrNum].sequence.seq + UnscoredMostHitedLocation - Indexer::Reference->chrs[chrNum].start_location, aSeq.length, aSeq.seq, aSeq.length, 100, L);
//	//cout << endl;
//	//while (true)
//	//{
//	//	cin >> index;for (unsigned i = index; i < index+100; ++i)
//	//	{
//	//		cout << Indexer::Reference->chrs[chrNum].sequence.seq[i - Indexer::Reference->chrs[chrNum].start_location];
//	//	}
//	//	cout << endl;
//	//}
//
//	return 0;
//}

	printf("Welcome to SNAPI version %s.\n\n", SNAPI_VERSION);
	try
	{
		if (argc < 2) {
			usage();
		}
		else if (strcmp(argv[1], "index") == 0) {
			buildIndex(argc - 2, argv + 2);
		}
		else if (strcmp(argv[1], "align") == 0) {
			runAlign(argc - 2, argv + 2);
		}
		else {
			fprintf(stderr, "Invalid command: %s\n\n", argv[1]);
			usage();
		}
	}
	catch (int err)
	{
		printError(err);
		return 1;
	}
	return 0;
}

static void iusage()
{
	fprintf(stderr,
		"Usage: snapi index <input.fa> <output-file> [<options>]\n"
		"Options:\n"
		"  -s  seed size (default: %u)\n"
		"  -h  hMax (default: %u)\n"
		"  -C  continue cutting (the <input.fa> will be seen as uncuted or unfinished cutted index)\n"
		"  -X  no cutting, build index the old fashioned way (won't cut the index after build)\n",
		seed::seedLen,
		Indexer::hMax);
	exit(1);
}

void buildIndex(int argc, const char **argv)
{
	if (argc < 2) {
		iusage();
	}

	const char *fastaFile = argv[0];
	const char *outputFile = argv[1];

	unsigned C = false, X = false;

	for (int n = 2; n < argc; n++) {
		if (strcmp(argv[n], "-s") == 0) {
			if (n + 1 < argc) {
				seed::seedLen = atoi(argv[n + 1]);
				n++;
			}
			else {
				iusage();
			}
		}
		else if (strcmp(argv[n], "-h") == 0) {
			if (n + 1 < argc) {
				Indexer::hMax = atof(argv[n + 1]);
				n++;
			}
			else {
				iusage();
			}
		}
		else if (strcmp(argv[n], "-C") == 0) {
			C = true;
		}
		else if (strcmp(argv[n], "-X") == 0) {
			X = true;
		}
		else {
			fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
			iusage();
		}
	}

	if (seed::seedLen < 16 || seed::seedLen > MAX_SEED_LENGTH) {
		// Right now there's some hard-coded stuff, like the seed warp table, that
		// does not work for seed lengths outside the range of 16-23.
		fprintf(stderr, "Seed length must be between 16 and %d\n",MAX_SEED_LENGTH);
		exit(1);
	}
	if (C == true && X == true)
	{
		fprintf(stderr, "Can't continue cutting when we don't even cut.");
		iusage();
	}
	Indexer TheIndexer;
	if (!C)
	{
		printf("Building index from FASTA file '%s'\n", fastaFile);
		printf("Building...\n");
		start = time(NULL);
		TheIndexer.readReference(fastaFile);
		TheIndexer.buildIndex();
		stop = time(NULL);
		printf("Hash table building succeeded! Time spent: %llds.\n Cutting Hash table...\n", stop - start);
		FILE* file;
		file = fopen(outputFile, "wb");
		if (file == nullptr)
		{
			throw - 102;
		}
		unsigned FinishedIndex = 0;
		fwrite(&FinishedIndex, sizeof(unsigned), 1, file);
		fclose(file);
	}
	else
	{
		printf("Reading index from index files '%s'\n", fastaFile);
		printf("Reading...\n");
		start = time(NULL);
		TheIndexer.loadFromFile(fastaFile);
		stop = time(NULL);
		printf("Hash table reading succeeded! Time spent: %llds.\n Cutting Hash table...\n", stop - start);
	}
	if (X)
	{
		TheIndexer.saveToFile(outputFile);
	}
	else
	{
		start = time(NULL);
		TheIndexer.cutIndex(Indexer::hMax, outputFile);
		stop = time(NULL);
		printf("Hash table cuted! Time spent: %llds.\n", stop - start);
	}
	//TheIndexer.cleanIndex();
}

static void ausage()
{
	fprintf(stderr,
		"Usage: snapi align <index-file> <inputFile> <outputFile(.sam)> [<options>]\n"
		"Options:\n"
		"- d   maximum edit distance allowed per read or pair(default: %u)\n"
		"- n   number of seeds to use per read(default: %u)\n"
		"- h   maximum hits to consider per seed(default: %u)\n"
		"- c   confidence threshold(default: %u)\n"
		"- t   thread number(default: 1)\n"
		,Aligner::dMax, Aligner::SeedsToTry, Indexer::hMax, Aligner::confidence
		);
	exit(1);
}

void runAlign(int argc, const char **argv)
{
	if (argc < 3) {
		ausage();
	}

	const char *indexFile = argv[0];
	const char *inputFile = argv[1];
	const char *outputFile = argv[2];

	for (int n = 3; n < argc; n++) {
		if (strcmp(argv[n], "-d") == 0) {
			if (n + 1 < argc) {
				Aligner::dMax = atoi(argv[n + 1]);
				if (Aligner::dMax > MAX_EDD)
				{
					fprintf(stderr,"Maximum edit distance allowed per read can be no bigger than %d.",MAX_EDD);
					exit(1);
				}
				n++;
			}
			else {
				ausage();
			}
		}
		else if (strcmp(argv[n], "-n") == 0) {
			if (n + 1 < argc) {
				Aligner::SeedsToTry = atof(argv[n + 1]);
				if (Aligner::SeedsToTry > MAX_SEED_NUM)
				{
					fprintf(stderr, "Number of seeds to use per read can be no bigger than %d.", MAX_SEED_NUM);
					exit(1);
				}
				n++;
			}
			else {
				ausage();
			}
		}
		else if (strcmp(argv[n], "-h") == 0) {
			if (n + 1 < argc) {
				Indexer::hMax = atof(argv[n + 1]);
				if (Indexer::hMax > MAX_HIT_NUM)
				{
					fprintf(stderr, "Maximum hits to consider per seed can be no bigger than %d.", MAX_HIT_NUM);
					exit(1);
				}
				n++;
			}
			else {
				ausage();
			}
		}
		else if (strcmp(argv[n], "-c") == 0) {
			if (n + 1 < argc) {
				Aligner::confidence = atof(argv[n + 1]);
				if (Aligner::confidence > MAX_EDD)
				{
					fprintf(stderr, "Confidence threshold can be no bigger than %d.", MAX_EDD);
					exit(1);
				}
				n++;
			}
			else {
				ausage();
			}
		}
		else if (strcmp(argv[n], "-t") == 0) {
			if (n + 1 < argc) {
				threadNum = atof(argv[n + 1]);
				if (threadNum > 8)
				{
					fprintf(stderr, "Thraed number can be no bigger than 8.");
					exit(1);
				}
				n++;
			}
			else {
				ausage();
			}
		}
		else {
			fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
			ausage();
		}
	}

	if (seed::seedLen < 16 || seed::seedLen > MAX_SEED_LENGTH) {
		// Right now there's some hard-coded stuff, like the seed warp table, that
		// does not work for seed lengths outside the range of 16-23.
		fprintf(stderr, "Seed length must be between 16 and %d\n", MAX_SEED_LENGTH);
		exit(1);
	}
	printf("Reading index...\n"); 
	start = time(NULL);
	Indexer TheIndexer;
	TheIndexer.loadFromFile(indexFile);
	stop = time(NULL);
	printf("Done reading index! Time spent: %llds.\nAligning...\n", stop - start); 
	start = time(NULL);
	Aligner TheAligner;
	TheAligner.singleAlign(inputFile, outputFile, threadNum);
	stop = time(NULL);
	unsigned speed = ((unsigned)TheAligner.SequenceAligned) / (unsigned)(stop - start);
	printf("All %u reads aligned! Time spent: %llds. Read speed: %u reads/s, hMax: %u.", TheAligner.SequenceAligned, stop - start, speed, Indexer::hMax);
	//TheIndexer.cleanIndex();
}
