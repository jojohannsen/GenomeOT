#include "suffix_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Usage()
{
	printf("Usage: chrcompare <file1> <start offset> <suffix tree string length> <file2> <segment size> <window size>\n");
	printf("\n");
	printf(" Reads in <suffix tree string length> characters from <file1> starting at <start offset>\n");
	printf(" then scans all of <file2> and in each <segment size> section, looks up each <window size>\n");
	printf(" string in the suffix tree, and counts how many times the string was found in each section.\n");
	printf(" \n");
	printf(" Outputs 1 line per section, comma delimited: <section offset>,<section count>,<count for each segment>\n");
}


// 
//  complement -- reverse complement a character
//
//    A => T
//    C => G
//    G => C
//    T => A
//
char complement( char cval )
{
	if (cval == 'A') return 'T';
	else if (cval == 'C') return 'G';
	else if (cval == 'G') return 'C';
	else if (cval == 'T') return 'A';
	return 'x';
}

void reverse_complement( char* buf, int window_size )
{
	int first_offset = 0;
	int last_offset = window_size - 1;
	char temp;
	// first_offset and last_offset move in toward center, then loop exits
	while (first_offset < last_offset)
	{
		// save the value we are about to overwrite, the "first_offset" value
		temp = *(buf + first_offset);
		// overwrite "first_offset" value with reverse complement of the "last_offset" value
		*(buf + first_offset) = complement(*(buf + last_offset));
		// overwrite "last_offset" value with reverse complement of our "temp" save of the original "first_offset" value
		*(buf + last_offset) = complement( temp );
		// move pointers in toward center
		first_offset++;
		last_offset--;
	}
	// in case where there is a lone character in the center, replace it with it's complement
	if (first_offset == last_offset)
	{
		*(buf + first_offset) = complement(*(buf + first_offset));
	}
}

int main(int argc, char* argv[])
{
	/* command line parameters */
	unsigned char* file1 = NULL;
	DBL_WORD start_offset = 0;
	DBL_WORD suffix_tree_string_length = 0;
	unsigned char* file2 = NULL;
	DBL_WORD segment_size = 0;
	DBL_WORD window_size = 0;

	/* internal data */
	SUFFIX_TREE* tree = NULL;;
	FILE* inFile1 = NULL;
	FILE* inFile2 = NULL;
	unsigned char* data_buffer = NULL;
	unsigned char* data_buffer2 = NULL;
	unsigned char* window = NULL;
	unsigned char* scanner = NULL;
	int offset = 0;
	int forward_count = 0;
	int backward_count = 0;
	int section_number = 0;
	int buckets_per_segment = 10;
	int* buckets;
	int* backward_buckets;
	int i = 0;
	DBL_WORD position;

	/* Set up parameters, validate */
	if (argc < 7) 
	{
		Usage();
		exit(0);
	}
	file1 = argv[1];
	start_offset = atol(argv[2]);
	suffix_tree_string_length = atol(argv[3]);
	file2 = argv[4];
	segment_size = atol(argv[5]);
	window_size = atol(argv[6]);
	buckets_per_segment = (int)(suffix_tree_string_length/segment_size);
	buckets = (int*)malloc(buckets_per_segment*sizeof(int));
	backward_buckets = (int*)malloc(buckets_per_segment*sizeof(int));

	/* open the file, scan to offset, read in characters, create suffix tree */
	inFile1 = fopen((const char*)file1, "r");
	if (inFile1 == NULL)
	{
		printf("File '%s' NOT FOUND.\n", file1);
		exit(0);
	}
	data_buffer = (unsigned char*)malloc(suffix_tree_string_length*sizeof(unsigned char));
	fseek( inFile1, start_offset, SEEK_CUR );
	fread( data_buffer, 1, suffix_tree_string_length, inFile1 );
	tree = ST_CreateTree((const char*)data_buffer, suffix_tree_string_length);

	/* open the second file, read in chunks, match each section against suffix tree */
	inFile2 = fopen((const char*)file2, "r");
	if (inFile2 == NULL)
	{
		printf("File '%s' NOT FOUND.\n", file2);
		exit(0);
	}
	data_buffer2 = (unsigned char*)malloc(segment_size);
	window = (unsigned char*)malloc( window_size + 1 );
	*(window + window_size) = 0;
	section_number = 0;

	while ( fread( data_buffer2, 1, segment_size, inFile2 ) == segment_size )
	{
		offset = 0;
		forward_count = 0;
		backward_count = 0;
		scanner = data_buffer2;
		for (i = 0; i < buckets_per_segment; i++) {
			*(buckets + i) = 0;
			*(backward_buckets + i) = 0;
		}
		while ( offset < segment_size )
		{
			strncpy( window, scanner, window_size );
			if ((position = ST_FindSubstring( tree, window, window_size )) != ST_ERROR)
			{
				forward_count++;
				i = (int)(position/segment_size);
				if (i < buckets_per_segment) {
					*(buckets + i) += 1;
				}
			}
			reverse_complement( window, window_size );
			if ((position = ST_FindSubstring( tree, window, window_size )) != ST_ERROR)
			{
				backward_count++;
				i = (int)(position/segment_size);
				if (i < buckets_per_segment) {
					*(backward_buckets + i) += 1;
				}
			}
			scanner += window_size;
			offset += window_size;
		}
		printf("%d,%d,%d", section_number++, forward_count, backward_count );
		for (i = 0; i < buckets_per_segment; i++) {
			printf(",%d", *(buckets + i));
		}
		for (i = 0; i < buckets_per_segment; i++) {
			printf(",%d", *(backward_buckets + i));
		}
		printf("\n");
		fflush(stdout);

	}
	free( data_buffer2 );
	free( data_buffer );
	return 0;
}
