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
	printf("Usage: centromere <file name> <window size> <overlap> [DAWG] [<min depth>-<max depth>] [<interval size>] \n");
	printf("\n");
	printf(" <window size> range %lu to %lu, and greater than 'overlap' value\n", MIN_WINDOW_SIZE, MAX_WINDOW_SIZE);
	printf(" <overlap> range %lu to %lu\n", MIN_OVERLAP, MAX_OVERLAP);
	printf(" [<interval size>] breaks depth range into chunks\n");
	printf(" [DAWG] removes nodes that have suffix links to nodes with same child counts.\n");
	printf(" [LEFT] removes nodes that are not left diverse.\n");
	printf("\n");
	printf("Breaks a file into overlapping windows, and for each window, prints the following values:");
	printf("\n");
	printf("The first 3 window values indicate the location of the start of the section:\n");
	printf("\n");
	printf("\tSource File Line Number\n");
	printf("\tSource File Line Number Offset\n");
	printf("\tSequence Offset\n");
	printf("\n");
	printf("The next values are pairs based on the optional parameters:\n");
	printf("\n");
	printf("\tNumber of Nodes in Suffix Tree\n");
	printf("\tNumber of Distinct Substrings in Suffix Tree\n");
	printf("\n");
	printf("Note: suffix trees allow many more calculations than this...\n");
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
 *    DBL_WORD* allocate_counts( generate_DAWG, min_depth, max_depth, interval_size )
 */
int number_counts = 0;
DBL_WORD* counts_memory = NULL;
DBL_WORD* counts_memory_scanner = NULL;
DBL_WORD file_line_number = 0;
DBL_WORD file_line_offset = 0;
DBL_WORD sequence_offset = 0;
/* duplicate the above 3, but since we backtrack, save these when the next chunk start location is first encountered,
   so the values are correct after the fseek back to this location.
 */
DBL_WORD NEXT_CHUNK_file_line_number = 0;
DBL_WORD NEXT_CHUNK_file_line_offset = 0;
DBL_WORD NEXT_CHUNK_sequence_offset = 0;

void counts_location_adjust()
{
	file_line_number = NEXT_CHUNK_file_line_number;
	file_line_offset = NEXT_CHUNK_file_line_offset;
	sequence_offset = NEXT_CHUNK_sequence_offset;
}

int counts_generate_DAWG = 0;
int counts_detect_left_diverse = 0;
DBL_WORD counts_min_depth = 0;
DBL_WORD counts_max_depth = 0;
DBL_WORD counts_interval_size = 0;

void print_counts()
{
	int i = 0;
	counts_memory_scanner = counts_memory;
	for (i = 0; i < number_counts; i++)
	{
		printf("%lu", *counts_memory_scanner++);
		if (i < (number_counts + 1))
		{
			printf(" ");
		}
	}
	printf("\n");
}

void allocate_counts( int generate_DAWG, int detect_left_diverse, DBL_WORD min_depth, DBL_WORD max_depth, DBL_WORD interval_size )
{
	counts_generate_DAWG = generate_DAWG;
	counts_detect_left_diverse = detect_left_diverse;
	counts_min_depth = min_depth;
	counts_max_depth = max_depth;
	counts_interval_size = interval_size;

	number_counts = 3;
	if (interval_size == 0)
	{
		number_counts += 2;
	}
	else
	{
		number_counts += (( max_depth - min_depth + interval_size - 1)/interval_size)*2;  /*  times 2 because we store 2 values for each interval */
	}
	counts_memory = (DBL_WORD*)malloc(number_counts*sizeof(DBL_WORD));
	memset(counts_memory, 0, number_counts*sizeof(DBL_WORD));
}


void generate_node_counts(SUFFIX_TREE* tree, NODE* node, 
						  long node_depth, long string_depth, 
						  long string_depth_start, long string_depth_end, 
						  DBL_WORD* node_count, DBL_WORD* substring_count, DBL_WORD* substring_millions_count)
{
   NODE* child_scanner = node->sons;
   long  start = node->edge_label_start, end;
   end     = get_node_label_end(tree, node);
   string_depth += (end - start + 1);

   if (node->ignore_NODE)
   {
	   return;
   }

   if ((string_depth_end == NO_DEPTH_LIMIT) || ((string_depth >= string_depth_start) && (string_depth <= string_depth_end)))
   {
	   *node_count += 1;
	   *substring_count += (end - start + 1);
	   if (*substring_count > 1000000)
	   {
		   *substring_millions_count += 1;
		   *substring_count -= 1000000;
	   }
   }
   while(child_scanner != NULL)
   {
      generate_node_counts(tree, child_scanner, node_depth + 1, string_depth, string_depth_start, string_depth_end, node_count, substring_count, substring_millions_count);
      child_scanner = child_scanner->right_sibling;
   }
}


void generate_left_diverse_nodes( NODE* node )
{
	NODE* child_scanner = node->sons;
	char left_char = 0;
	int is_left_diverse = 1;
	if (child_scanner == NULL) 
	{
		node->is_left_diverse = 1;
	}
	else
	{	
		while (child_scanner != NULL)
		{
			generate_left_diverse_nodes( child_scanner );
			if (left_char == 0)
			{
				left_char = child_scanner->left_char;
			} 
			else if (left_char != child_scanner->left_char)
			{
				is_left_diverse = 0;
			}
			if (child_scanner->is_left_diverse == 0)
			{
				is_left_diverse = 0;
			}
			child_scanner = child_scanner->right_sibling;
		}
		node->is_left_diverse = is_left_diverse;
	}
	if (!node->is_left_diverse)
	{
		node->ignore_NODE = 1;
	}
}

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

void generate_DAWG_nodes( NODE* node )
{
	NODE* child_scanner = node->sons;
	/* If this node has a suffix link to another NODE with the same leaf count, set the "ignore_NODE" flag */
	if (( node->suffix_link != NULL ) && ( node->leaf_count == node->suffix_link->leaf_count ))
	{
		node->ignore_NODE = 1;
	}
	else
	{
		while (child_scanner != NULL)
		{
			generate_DAWG_nodes( child_scanner );
			child_scanner = child_scanner->right_sibling;
		}
	}
}

void generate_counts( SUFFIX_TREE* tree )
{
	/* counts that accumulate during traversal */
	DBL_WORD node_count = 0;
	DBL_WORD substring_count = 0;
	DBL_WORD substring_millions_count = 0;
	DBL_WORD interval_start_depth = 0;
	DBL_WORD interval_end_depth = 0;

	if (counts_generate_DAWG)
	{
		generate_leaf_counts( tree->root );
		generate_DAWG_nodes( tree->root );
	}
	if (counts_detect_left_diverse)
	{
		generate_left_diverse_nodes( tree->root );
	}
	if ((counts_max_depth == NO_DEPTH_LIMIT) || (counts_interval_size == 0))
	{
		generate_node_counts( tree, tree->root, 0, 0, counts_min_depth, counts_max_depth, &node_count, &substring_count, &substring_millions_count );
		*counts_memory_scanner++ = node_count;
		*counts_memory_scanner++ = substring_count;
	}
	else
	{
		interval_start_depth = counts_min_depth;
		interval_end_depth = counts_min_depth + counts_interval_size - 1;
		while (interval_start_depth < counts_max_depth)
		{
			node_count = substring_count = 0;
			generate_node_counts( tree, tree->root, 0, 0, interval_start_depth, interval_end_depth, &node_count, &substring_count, &substring_millions_count );
			*counts_memory_scanner++ = node_count;
			*counts_memory_scanner++ = substring_millions_count;
			interval_start_depth += counts_interval_size;
			interval_end_depth += counts_interval_size;
			interval_end_depth = (interval_end_depth > counts_max_depth) ? counts_max_depth : interval_end_depth;
		}
	}
}

/*
 *  Updates:
 *
 *    file_line_number -- any time '\n' is encountered, this is incremented
 *    file_line_offset -- resets each time '\n' encountered, ignores '\r' in line_offset
 *    sequence_offset -- increments each time acgt or ACGT encountered, converts acgt to ACGT when put in buffer
 *
 *  Internally, has in_comment, set to true when '>' encountered, resets with '\n' (all characters in between ignored)
 */
DBL_WORD counts_fread( char* data_buffer, DBL_WORD number_bytes, FILE* file, DBL_WORD overlap )
{
	DBL_WORD bytes_read = 0;
	int in_comment = 0;
	int cval = 0;
	int file_offset = 0;

	/* save chunk values for when we need to print them along with chunk stats */
	memset(counts_memory, 0, number_counts*sizeof(DBL_WORD));
	counts_memory_scanner = counts_memory;
	*counts_memory_scanner++ = file_line_number;
	*counts_memory_scanner++ = file_line_offset;
	*counts_memory_scanner++ = sequence_offset;

	while ((bytes_read < number_bytes) || (in_comment))
	{
		cval = fgetc( file );

		if (feof(file)) break;
		if (cval == '\n')
		{
			in_comment = 0;
			file_line_number++;
			file_line_offset = 0;
		}
		else if (cval == '>')
		{
			in_comment = 1;
		}
		if (!in_comment)
		{
			if (cval != '\r')
			{
				file_line_offset++;
			}
			if (cval == 'a') cval = 'A';
			else if (cval == 'c') cval = 'C';
			else if (cval == 'g') cval = 'G';
			else if (cval == 't') cval = 'T';
			if (cval == 'N') 
			{
				sequence_offset++;
			}
			if ((cval == 'A') || (cval == 'C') || (cval == 'G') || (cval == 'T'))
			{
				*data_buffer++ = (char)cval;
				bytes_read++;
				sequence_offset++;
				if (bytes_read == (number_bytes - overlap))
				{
					NEXT_CHUNK_file_line_number = file_line_number;
					NEXT_CHUNK_file_line_offset = file_line_offset;
					NEXT_CHUNK_sequence_offset = sequence_offset;
				}

			}
		}
	}
	return bytes_read;
}

void print_counts_header( int generate_DAWG, DBL_WORD min_depth, DBL_WORD max_depth, DBL_WORD interval_size )
{
	int first_interval = 0;
	int last_interval = 0;

	printf("# LineNo, LineOffset, SeqOffset, ");
	if (interval_size == 0)
	{
		printf("Nodes, Substrings, ");
	}
	else
	{
		first_interval = min_depth;
		last_interval = min_depth + interval_size - 1;
		while ( first_interval <= max_depth )
		{
			printf("Nodes(%d,%d), Substrings(%d,%d), ", first_interval, last_interval, first_interval, last_interval);
			first_interval += interval_size;
			last_interval += interval_size;
			if (last_interval > max_depth)
			{
				last_interval = max_depth;
			}
		}
	}
	if (generate_DAWG)
	{
		printf("DAWG_Nodes, DAWG_Substrings, ");
	}
	printf("\n");
}

/* command line parameter support methods */
#define valid_range( value, mnv, mxv ) ((value>=mnv)&&(value<=mxv))

void extract_range( DBL_WORD* range_min, DBL_WORD* range_max, int argc, char* argv[] )
{
	int i = 0;
	char* range_str = NULL;
	char tokenizer_buffer[100];
	char* token;

	for (i = 0; i < argc; i++)
	{
		range_str = strstr( argv[i], "-" );
		if (range_str != NULL)
		{
			strcpy(tokenizer_buffer, argv[i]);
			token = strtok(tokenizer_buffer,"-");
			*range_min = atol(token);
			token = strtok(NULL, "-");
			*range_max = atol(token);
			return;
		}
	}
}

void extract_flag( int* flag, const char* str, int argc, char* argv[] )
{
	int i = 0;

	for (i = 0; i < argc; i++) 
	{
		if (strcmp( argv[i], str ) == 0)
		{
			*flag = 1;
			return;
		}
	}
}

void extract_long( DBL_WORD* long_val, int min_offset_to_check, int argc, char* argv[] )
{
	int i = 0;
	DBL_WORD test_val = 0;

	for (i = min_offset_to_check; i < argc; i++)
	{
		if (strstr( argv[i], "-" ) == NULL)
		{
			test_val = atol( argv[i] );
			if (test_val > 0)
			{
				*long_val = test_val;
				return;
			}
		}
	}
}

int main(int argc, char* argv[])
{
	/* command line parameters */
	unsigned char* file_name = NULL;
	DBL_WORD window_size = 0;
	DBL_WORD overlap = 0;
	int generate_DAWG = 0;
	int detect_left_diverse = 0;
	DBL_WORD min_depth = NO_DEPTH_LIMIT;
	DBL_WORD max_depth = NO_DEPTH_LIMIT;
	DBL_WORD interval_size = 0;

	/* internal data */
	SUFFIX_TREE* tree = NULL;;
	FILE* file = NULL;
	unsigned char* data_buffer = NULL;
	DBL_WORD* counts = NULL;

	/* Set up parameters, validate */
	if (argc < 4) 
	{
		Usage();
		exit(0);
	}
	file_name = argv[1];
	window_size = atol(argv[2]);
	overlap = atol(argv[3]);

	if (!valid_range( window_size, MIN_WINDOW_SIZE, MAX_WINDOW_SIZE ) || !valid_range( overlap, MIN_OVERLAP, MAX_OVERLAP ) || (window_size < overlap))
	{
		Usage();
		exit(0);
	}

	extract_range( &min_depth, &max_depth, argc, argv );
	extract_flag( &generate_DAWG, "DAWG", argc, argv );
	extract_flag( &detect_left_diverse, "LEFT", argc, argv );

	/* 
	 * Parameters:  centromere <file name> <window size> <overlap> [DAWG] [<min depth>-max depth>] [<interval size>]
	 * Offsets:     0          1           2             3         >3     >3                       >3
	 */
	int min_offset_to_check = 4;
	extract_long( &interval_size, min_offset_to_check, argc, argv );

	/* everything checks out, run the algorithm */

	/* open the file, read in chunks of 'window_size', create suffix tree, generate counts and print them, then back track in file by 'overlap', then repeat */
	file = fopen((const char*)file_name, "r");
	if (file == NULL)
	{
		printf("File '%s' NOT FOUND.\n", file_name);
		exit(0);
	}
	data_buffer = (unsigned char*)malloc(window_size*sizeof(unsigned char));
	allocate_counts( generate_DAWG, detect_left_diverse, min_depth, max_depth, interval_size );

	/* read in chunks of 'window_size', create suffix tree, generate counts, print them, back track by 'overlap' */
	print_counts_header( generate_DAWG, min_depth, max_depth, interval_size );
	while (counts_fread( data_buffer, window_size, file, overlap ) == window_size)
	{
		tree = ST_CreateTree((const char*)data_buffer, window_size);
		generate_counts( tree );
		print_counts();
		fseek( file, -overlap, SEEK_CUR );
		ST_DeleteTree( tree );
		counts_location_adjust( overlap );
	}
	free( data_buffer );
	free( counts );
	return 0;
}
