#include "suffix_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* constants related to command line parameter values */
#define NO_DEPTH_LIMIT -1
#define MIN_WINDOW_SIZE 10L
#define MAX_WINDOW_SIZE 100000000L
#define MIN_OVERLAP 0L
#define MAX_OVERLAP 1000000L

void Usage()
{
	printf("Usage: st_scan <suffix tree file name> <file to scan> <scan size>\n");
	printf("\n");
	printf(" <scan size> is a fixed window size to check against suffix tree\n");
}

char rc( char cval )
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
	while (first_offset < last_offset)
	{
		temp = *(buf + first_offset);
		*(buf + first_offset) = rc(*(buf + last_offset));
		*(buf + last_offset) = rc( temp );
		first_offset++;
		last_offset--;
	}
	if (first_offset == last_offset)
	{
		*(buf + first_offset) = rc(*(buf + first_offset));
	}
}


/* counts class methods:
 *
 * All functions include parameters which include the number of values to extract from suffix tree: 
 *
 *    3 location values (globals), plus
 *
 *    generate_DAWG => 2 values
 *    min_depth, max_depth, interval_size => if interval_size=0, 2 values, othersize ((max_depth - min_depth + interval_size - 1)/interval_size)*2 value
 *        
 */
int number_counts = 0;
DBL_WORD file_line_number = 0;
DBL_WORD file_line_offset = 0;
DBL_WORD sequence_offset = 0;
/* duplicate the above 3, but since we backtrack, save these when the next chunk start location is first encountered,
   so the values are correct after the fseek back to this location.
 */
DBL_WORD NEXT_CHUNK_file_line_number = 0;
DBL_WORD NEXT_CHUNK_file_line_offset = 0;
DBL_WORD NEXT_CHUNK_sequence_offset = 0;


void generate_leaf_counts( NODE* node )
{
	/* For this, I have to add a "leaf_count" column to NODE */
	int lc = 0;
	NODE* child_scanner = node->sons;

	if (child_scanner == NULL) 
	{
		lc = 1;
	}
	else
	{
		while (child_scanner != NULL)
		{
			generate_leaf_counts( child_scanner );
			lc += child_scanner->leaf_count;
			child_scanner = child_scanner->right_sibling;
		}
	}
	node->leaf_count = lc;
}




int main(int argc, char* argv[])
{
	/* command line parameters */
	unsigned char* st_file_name = NULL;
	DBL_WORD st_file_size = 0;
	unsigned char* scan_file_name = NULL;
	char* scan_buffer = NULL;
	DBL_WORD window_size = 0;
	DBL_WORD position = 0;
	int foundCount = 0;
	int notFoundCount = 0;
	int found = 0;
	int forwardCount = 0;
	int backwardCount = 0;

	/* internal data */
	SUFFIX_TREE* tree = NULL;
	FILE* file = NULL;
	FILE* fileToScan = NULL;
	unsigned char* data_buffer = NULL;

	/* Set up parameters, validate */
	if (argc != 4) 
	{
		Usage();
		exit(0);
	}
	st_file_name = argv[1];
	scan_file_name = argv[2];
	window_size = atol(argv[3]);


	file = fopen((const char*)st_file_name, "r");
	if (file == NULL)
	{
		printf("File '%s' NOT FOUND.\n", st_file_name);
		exit(0);
	}
	fseek(file, 0L, SEEK_END);
	st_file_size = ftell( file );
	fseek( file, 0L, SEEK_SET );

	data_buffer = (unsigned char*)malloc(st_file_size*sizeof(unsigned char));
	fread( data_buffer, st_file_size, 1, file );
        fclose( file );
	tree = ST_CreateTree((const char*)data_buffer, st_file_size);

	/* from here, read sections of the file to scan, count the number of forward and
         * reverse matches.
         *
         * when scanning is done, print the results as "forward,backward"
         */
 
        fileToScan = fopen(scan_file_name, "r");
	scan_buffer = (char*)malloc( window_size + 1 );
	while (fread( scan_buffer, window_size, 1, fileToScan ) == 1) 
	{
		found = 0;

		if ((position = ST_FindSubstring( tree, scan_buffer, window_size )) != ST_ERROR)
		{
			forwardCount++;
			found = 1;
		}
		reverse_complement( scan_buffer, window_size );

		if ((position = ST_FindSubstring( tree, scan_buffer, window_size )) != ST_ERROR)
		{
			backwardCount++;
			found = 1;
		}
		if (found) 
		{
			foundCount++;
		}
		else
		{
			notFoundCount++;
		}
	}
	fclose( fileToScan );
	printf("%s,%s,%d,%d,%d,%d\n", st_file_name, scan_file_name, forwardCount, backwardCount, foundCount, notFoundCount);

	free( data_buffer );
	return 0;
}
