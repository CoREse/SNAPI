#include <iostream>
#include "seed.h"
#include "index.h"
#include "sequence.h"
#include <stdio.h>
#include <time.h>
#include <thread>
#include <vector>
#include <mutex>
#include "base.h"
#include "NASeq.h"

using namespace std;


const char *SNAP_VERSION = "0.11";
time_t start, stop;

unsigned threadNum=1;

vector<thread> threads;
mutex mtxR,mtxW;

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

index MyIndex;

void buildIndex(int argc, const char **argv);
void runAlign(int, const char **);

int main(int argc, const char **argv)
{
	FILE *ifile = fopen("chr21.fa", "r");
	if (ifile == nullptr)
	{
		throw - 102;//file pointer is null
		return false;
	}
	NASeq aseq(ifile);
	return 0;
	cout << sizeof(aseq);
	printf("Welcome to SNAP version %s.\n\n", SNAP_VERSION);
	cin >> start;
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
	///*MyIndex.build("chr21.fa");
	//cout << MyIndex.Index[1].key << ' ' << MyIndex.Index[1].value << endl;
	//cout << (MyIndex.Index + 13)->key << endl;
	//seed aseed("aaaaaaaaaaaaaaaaaaaa");
	//unsigned * a;
	//unsigned num=MyIndex.lookup(aseed.toNum(),a);
	//cout << num;
	//MyIndex.save2file("chr21.index");
	//free(a);*/
	//MyIndex.readFromFile("chr21.index");
	////cout << MyIndex.getSize();
	//sequence seq;
	//FILE* fq = fopen("chr21.fq", "r");
	//if (nullptr == fq)
	//{
	//	cout << "error";
	//	return 0;
	//}
	//while(seq.readNext(fq))
	//seq.align(MyIndex);
	//if (sequence::outputFile != nullptr)
	//{
	//	fclose(sequence::outputFile);
	//}
	//fclose(fq);
	return 0;
}

static void iusage()
{
	fprintf(stderr,
		"Usage: snap index <input.fa> <output-file> [<options>]\n"
		"Options:\n"
		"  -s  seed size (default: %d)\n"
		"  -h  hash table slack (default: %.1f)\n",
		seed::length,
		index::slack);
	exit(1);
}

void buildIndex(int argc, const char **argv)
{
	if (argc < 2) {
		iusage();
	}

	const char *fastaFile = argv[0];
	const char *outputFile = argv[1];

	for (int n = 2; n < argc; n++) {
		if (strcmp(argv[n], "-s") == 0) {
			if (n + 1 < argc) {
				seed::length = atoi(argv[n + 1]);
				n++;
			}
			else {
				iusage();
			}
		}
		else if (strcmp(argv[n], "-h") == 0) {
			if (n + 1 < argc) {
				index::slack = atof(argv[n + 1]);
				n++;
			}
			else {
				iusage();
			}
		}
		else {
			fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
			iusage();
		}
	}

	if (seed::length < 16 || seed::length > MAX_SEED_LENGTH) {
		// Right now there's some hard-coded stuff, like the seed warp table, that
		// does not work for seed lengths outside the range of 16-23.
		fprintf(stderr, "Seed length must be between 16 and %d\n",MAX_SEED_LENGTH);
		exit(1);
	}

	printf("Building index from FASTA file '%s'\n", fastaFile);
	printf("Hash table slack %lf\nBuilding...\n", index::slack);
	start = time(NULL);
	MyIndex.build(fastaFile);
	stop = time(NULL);
	printf("Hash table building succeeded! Time spent: %llds.\n Saving data into files...\n", stop - start); 
	start = time(NULL);
	MyIndex.save2file(outputFile);
	stop = time(NULL);
	printf("Data saved! Time spent: %llds.\n", stop - start); 
}

static void ausage()
{
	fprintf(stderr,
		"Usage: snap align <index-file> <inputFile> <outputFile(.sam)> [<options>]\n"
		"Options:\n"
		"- d   maximum edit distance allowed per read or pair(default: 8)\n"
		"- n   number of seeds to use per read(default: 25)\n"
		"- h   maximum hits to consider per seed(default: 250)\n"
		"- c   confidence threshold(default: 2)\n"
		"- t   thread number(default: 1)\n"
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
				sequence::dMax = atoi(argv[n + 1]);
				if (sequence::dMax > MAX_EDD)
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
				sequence::SeedsToTry = atof(argv[n + 1]);
				if (sequence::SeedsToTry > MAX_SEED_NUM)
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
				sequence::hMax = atof(argv[n + 1]);
				if (sequence::hMax > MAX_HIT_NUM)
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
				sequence::c = atof(argv[n + 1]);
				if (sequence::c > MAX_EDD)
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
					fprintf(stderr, "Thraed number can be no bigger than %d.", MAX_EDD);
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

	if (seed::length < 16 || seed::length > MAX_SEED_LENGTH) {
		// Right now there's some hard-coded stuff, like the seed warp table, that
		// does not work for seed lengths outside the range of 16-23.
		fprintf(stderr, "Seed length must be between 16 and %d\n", MAX_SEED_LENGTH);
		exit(1);
	}
	printf("Reading index...\n"); 
	start = time(NULL);
	MyIndex.readFromFile(indexFile);
	stop = time(NULL);
	for (int i = 0; i != 256; ++i)
	{
		sequence::samName[i] = outputFile[i];
		if (outputFile[i] == '\0') break;
	}
	printf("Done reading index! Time spent: %llds.\nAligning...\n", stop - start); 
	FILE* fq = fopen(inputFile, "r");
	if (nullptr == fq)
	{
		cout << "Read File Error!"<<endl;
		exit(0);
	}
	start = time(NULL);
	for (int i = 0; i != threadNum; ++i)
	{
		threads.push_back(thread([fq](){
			sequence read, *bestS=nullptr;
			bool suc;
			do
			{
				mtxR.lock();
				suc = read.readNext(fq);
				mtxR.unlock();
				if (suc)
				{
					read.align(MyIndex, bestS);
					mtxW.lock();
					bestS->writeResult(MyIndex);
					mtxW.unlock();
					/*if (sequence::SequenceAligned % 10000 == 0)
					{
						printf("%u reads aligned\n", sequence::SequenceAligned);
					}*/
				}
			} while (suc);
		}));
	}
	for (auto& athread : threads){
		athread.join();
	}
	/*sequence read, *bestS = nullptr;
	bool suc;
	do
	{
		suc = read.readNext(fq);
		if (suc)
		{
			read.align(MyIndex, bestS);
			bestS->writeResult(MyIndex);
			if (sequence::SequenceAligned % 10000 == 0)
			{
				printf("%u reads aligned\n", sequence::SequenceAligned);
			}
		}
	} while (suc);*/
	stop = time(NULL);
	long speed = ((long)sequence::SequenceAligned) / (long)(stop - start);
	printf("All %u reads aligned! Time spent: %llds. Read speed: %ld reads/s.", sequence::SequenceAligned, stop - start, speed); 
	fclose(fq);
	if (sequence::outputFile != nullptr)
	{
		fclose(sequence::outputFile);
	}
}