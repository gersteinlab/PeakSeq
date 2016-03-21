#include <stdio.h>
#include <stdlib.h>
#include "mapped_read_tools.h"
#include <vector>
#include <ctype.h>
#include "signal_track_tools.h"
#include "annot_region_tools.h"
#include "genomics_coords.h"
#include "utils.h"
#include "nomenclature.h"
#include "rng.h"
#include "seed_manager.h"
#include <string.h>
#include <algorithm>
#include "ansi_string.h"
#include "annot_region_tools.h"
#include "mapped_read_tools.h"

using namespace std; 

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

#define MAX(x, y) ((x)>(y))?(x):(y)
#define MIN(x, y) ((x)<(y))?(x):(y)

vector<char*>* sort_bucket_read_lines(char* bucket_fp)
{
	// Load the reads.
	vector<char*>* bucket_read_lines = buffer_file(bucket_fp);	
	vector<int>* read_starts = new vector<int>();
	vector<t_read_line_sorting_info*>* sorting_info_list = new vector<t_read_line_sorting_info*>();
	for(int i_read = 0; i_read < (int)bucket_read_lines->size(); i_read++)
	{
		int cur_read_start = 0;
		sscanf(bucket_read_lines->at(i_read), "%*s %*s %d", &cur_read_start);

		t_read_line_sorting_info* cur_line_info = new t_read_line_sorting_info();
		cur_line_info->start = cur_read_start;
		cur_line_info->read_line = bucket_read_lines->at(i_read);

		sorting_info_list->push_back(cur_line_info);
	} // i_read loop.

	sort(sorting_info_list->begin(), sorting_info_list->end(), sort_read_line_info);
	vector<char*>* sorted_bucket_read_lines = new vector<char*>();

	for(int i_read = 0; i_read < (int)sorting_info_list->size(); i_read++)
	{
		sorted_bucket_read_lines->push_back(sorting_info_list->at(i_read)->read_line);

		delete sorting_info_list->at(i_read);
	} // i_read loop.

	delete(sorting_info_list);
	delete(read_starts);

	delete(bucket_read_lines);

	return(sorted_bucket_read_lines);
}

// Following is for sorting the mapped reads offline.
bool sort_read_line_info(t_read_line_sorting_info* info1, t_read_line_sorting_info* info2)
{
	return(info1->start < info2->start);
}

void load_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_sam = NULL;
	if(strcmp(sam_fp, "stdin") == 0)
	{
		f_sam = stdin;
	}
	else
	{
		f_sam = open_f(sam_fp, "r");
	}

	char cur_seq_id[1000];
	char cur_seq_nucs[1000];
	char phred_quality_str[1000];
	while(1)
	{
		char* cur_sam_line = getline(f_sam);
		if(cur_sam_line == NULL)
		{
			break;
		}

		// Parse the sam read line.
		if(sscanf(cur_sam_line, "%s %*d %*s %*d %*d %*s %*s %*d %*d %s %s", cur_seq_id, cur_seq_nucs, phred_quality_str) != 3)
		{
			fprintf(stderr, "Could not parse %s\n", cur_sam_line);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = t_string::copy_me_str(cur_seq_id);
		new_sequenced_read->nucs = t_string::copy_me_str(cur_seq_nucs);
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = t_string::copy_me_str(phred_quality_str);
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	// If the sam file is not stdin, close it.
	if(strcmp(sam_fp, "stdin") != 0)
	{
		fclose(f_sam);
	}
}

// FASTQ LOADING INTERFACE:
// Following are for fastq loading from mapped read files.
void load_sequenced_reads_per_fastq(char* fastq_fp, vector<t_sequenced_read*>* sequenced_reads)
{
	// Use the piped input.
	FILE* f_fastq = NULL;
	if(strcmp(fastq_fp, "stdin") == 0)
	{
		f_fastq = stdin;
	}
	else
	{
		f_fastq = open_f(fastq_fp, "r");
	}

	while(1)
	{
		char* cur_seq_id = getline(f_fastq);
		if(cur_seq_id == NULL)
		{
			break;
		}

		char* cur_seq_nucs = getline(f_fastq);
		if(cur_seq_nucs == NULL)
		{
			fprintf(stderr, "Could not read the sequence for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_id = getline(f_fastq);
		if(cur_seq_qual_id == NULL)
		{
			fprintf(stderr, "Could not read the quality id for %s.\n", cur_seq_id);
			exit(0);
		}

		char* cur_seq_qual_str = getline(f_fastq);
		if(cur_seq_qual_str == NULL)
		{
			fprintf(stderr, "Could not read the quality str for %s.\n", cur_seq_id);
			exit(0);
		}

		// Add the new sequenced read.
		t_sequenced_read* new_sequenced_read = new t_sequenced_read();
		new_sequenced_read->id = cur_seq_id;
		new_sequenced_read->nucs = cur_seq_nucs;
		new_sequenced_read->mapping_info = NULL;
		new_sequenced_read->quality_str = cur_seq_qual_str;
		sequenced_reads->push_back(new_sequenced_read);

		if(sequenced_reads->size() % 100000 == 0)
		{
			fprintf(stderr, "Processing %ld. read.          \r", sequenced_reads->size());
		}
	} // file reading loop.

	if(strcmp(fastq_fp, "stdin") != 0)
	{
		fclose(f_fastq);
	}
}
void delete_sequenced_reads(vector<t_sequenced_read*>* sequenced_reads)
{
	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		delete [] sequenced_reads->at(i_r)->id;
		delete [] sequenced_reads->at(i_r)->nucs;
		delete [] sequenced_reads->at(i_r)->quality_str;

		delete sequenced_reads->at(i_r);
	} // i_r loop.

	delete sequenced_reads;
}

void dump_fastq(vector<t_sequenced_read*>* sequenced_reads, char* op_fastq_fp)
{
	FILE* f_op = open_f(op_fastq_fp, "w");

	for(int i_r = 0; i_r < (int)sequenced_reads->size(); i_r++)
	{
		fprintf(f_op, "@%s\n%s\n+%s\n%s\n", 
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->nucs,
			&(sequenced_reads->at(i_r)->id[1]), 
			sequenced_reads->at(i_r)->quality_str);
	} // i_r loop.

	fclose(f_op);
}

void load_mapped_sequenced_reads_per_SAM(char* sam_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_ELAND(char* eland_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_tagAlign(char* tagalign_fp, vector<t_sequenced_read*>* sequenced_reads);
void load_mapped_sequenced_reads_per_bowtie(char* bowtie_fp, vector<t_sequenced_read*>* sequenced_reads);
// ABOVE IS THE FASTQ LOADING INTERFACE:

// Generic preprocessing function for mapped read files.
void preprocess_mapped_reads_file(char* mrf_fp, char* parsed_reads_op_dir, void (preprocess_mapped_read_line)(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str),
	bool dump_read_id)
{
    // Divide SAM output with respect to chromosomes.
    FILE* f_mrf = NULL;
	t_file_buffer* mrf_file_buffer = NULL;
	if(strcmp(mrf_fp, "stdin") == 0)
	{
		f_mrf = stdin;
	}
	else
	{
		mrf_file_buffer = load_file(mrf_fp);
	}

	if(f_mrf == NULL && mrf_file_buffer == NULL)
	{
		fprintf(stderr, "mapped read file pointer and file buffer are both NULL for %s\n", mrf_fp);
		return;
	}

    //char cur_line[100000];
    int n_frags = 0;
    //int n_total_frags = 0;

    vector<FILE*>* frag_f_ptrs = new vector<FILE*>();

	// Check chromosome id's list file.
	char chr_ids_fp[100000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", parsed_reads_op_dir);

	vector<char*>* chr_ids = NULL;
	if(check_file(chr_ids_fp))
	{
		chr_ids = buffer_file(chr_ids_fp);

		fprintf(stderr, "Found chromosome id's @ %s, pooling.\n", chr_ids_fp);

		// Open the files for appending.
		for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
		{
			char new_fn[1000];
			sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
			fprintf(stderr, "Opening %s for pooling.\n", new_fn);
			if(!check_file(new_fn))
			{
				fprintf(stderr, "Could not open %s\n", new_fn);
				open_f(chr_ids_fp, "w");
				exit(0);
			}
			else
			{
				frag_f_ptrs->push_back(open_f(new_fn, "a"));
			}
		} // i_chr loop.
	}
	else
	{
		// The chromosomes will be added now.
		chr_ids = new vector<char*>();
	}

	while(1)
	{
		//char* cur_line = getline(f_mrf);
		char* cur_line = NULL;
		if(mrf_file_buffer != NULL)
		{
			cur_line = getline_per_file_buffer(mrf_file_buffer);
		}
		else if(f_mrf != NULL)
		{
			cur_line = getline(f_mrf);
		}

		if(cur_line == NULL)
		{
			break;
		}
		
		// Load the mapping info based on the file type.
		char chrom[1000];
		char read_id[1000];
		int chr_index;
		int sequenced_length;
		char strand_char; 
		char mapping_quality_str[20000];
		preprocess_mapped_read_line(cur_line,
									read_id,
									chrom, 
									chr_index, sequenced_length, 
									strand_char, 
									mapping_quality_str);

		// Make sure that the line is valid.
		if(chr_index >= 1 &&
			chrom[0] != 0)
		{
			// Normalize the chromosome id to comply with the 
			normalize_chr_id(chrom);

			// Get the chromosome index.
			int i_chr = t_string::get_i_str(chr_ids, chrom);

			// If the chromosome does not exist in the list opened and accumulated so far, add the id to the list and also open the processed read file.
			if(i_chr == (int)chr_ids->size())
			{
				// Add the chromosome id.
				chr_ids->push_back(t_string::copy_me_str(chrom));
				i_chr = t_string::get_i_str(chr_ids, chrom);

				char new_fn[10000];
				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chrom);

				// Does the file exist? If so, use the file, do not overwrite.
				frag_f_ptrs->push_back(open_f(new_fn, "w"));

				fprintf(stderr, "Added %s\n", chrom);
			}

			FILE* cur_frag_file = frag_f_ptrs->at(i_chr);

			if(cur_frag_file == NULL)
			{
					//printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
			}
			else
			{
				if(dump_read_id)
				{
					fprintf(cur_frag_file, "%s %s %c %d\n", read_id, mapping_quality_str, strand_char, chr_index);
				}
				else
				{
					fprintf(cur_frag_file, "%s %c %d\n", mapping_quality_str, strand_char, chr_index);
				}
				n_frags++;
			}
		} // check if the line corresponds to a valid mapped nucleotide.

		delete [] cur_line;
	} // file reading loop.

	// (Re)Dump the chromosome id list.
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for(int i_chr = 0; i_chr< (int)chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.

	fclose(f_chrs);

	// Close fragment file pointers.
	for(int i_f = 0; i_f < (int)frag_f_ptrs->size(); i_f++)
	{
		fclose(frag_f_ptrs->at(i_f));
	}

	// Unload/close the mapped read file.
	if(f_mrf != NULL)
	{
		fclose(f_mrf);
	}
	else if(mrf_file_buffer != NULL)
	{
		unload_file(mrf_file_buffer);
	}
}

void preprocess_tagAlign_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	int chr_start_index;
	int chr_end_index;
	char strand_sign;

	if(sscanf(cur_line, "%s %d %d %*s %*d %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	{
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		//chr_start_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		//chr_end_index += (CODEBASE_START_BASE - tagAlign_START_BASE);
		chr_start_index = translate_coord(chr_start_index, TAGALIGN_COORDS::start_base, CODEBASE_COORDS::start_base);
		chr_end_index = translate_coord(chr_end_index, TAGALIGN_COORDS::end_base, CODEBASE_COORDS::end_base);

		// Set quality to all matches.
		sprintf(cigar_str, "%dM", chr_end_index-chr_start_index+1);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = chr_end_index-chr_start_index+1;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_preprocessed_LH_GFF3_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int _chr_index;
	char strand;

	//if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %s", &flag, chrom, &_chr_index, mapping_quality_str, fragment, phred_quality_str) == 6)
	// X       14705460        35M     +
	if(sscanf(cur_line, "%s %d %s %c", chrom, &_chr_index, mapping_quality_str, &strand) == 4)
	{
		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		chr_index = _chr_index;

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand == '-')
		{
			strand_char = 'R';
		}

			chr_index = _chr_index;
			sequenced_length = 0;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_SAM_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* cigar_str)
{
	// Skip the comment and headers.
	if(cur_line[0] == '@')
	{
		chrom[0] = 0;
		chr_index = 0;
		return;
	}

	int flag;
	char flag_str[100];
	int _chr_index;
	char _chr_index_str[100];
	char fragment[100000];
	char phred_quality_str[100000];

	//t_string_tokens* cur_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	//if(sscanf(cur_line, "%s %d %s %d %*s %s %*s %*s %*s %s %s", read_id, &flag, chrom, &_chr_index, cigar_str, fragment, phred_quality_str) == 7)
	//if(cur_tokens->size() >= 11)
	if(sscanf(cur_line, "%[^'\t'] %[^'\t'] %[^'\t'] %[^'\t'] %*[^'\t'] %[^'\t'] %*[^'\t'] %*[^'\t'] %*[^'\t'] %[^'\t'] %[^'\t']", read_id, flag_str, chrom, _chr_index_str, cigar_str, fragment, phred_quality_str) == 7)
	{
		//t_string::copy(read_id, cur_tokens->at(0)->str());
		//flag = atoi(cur_tokens->at(1)->str());
		//t_string::copy(chrom, cur_tokens->at(2)->str());
		//_chr_index = atoi(cur_tokens->at(3)->str());
		//t_string::copy(cigar_str, cur_tokens->at(5)->str());
		//t_string::copy(fragment, cur_tokens->at(9)->str());
		//t_string::copy(phred_quality_str, cur_tokens->at(10)->str());

		_chr_index = atoi(_chr_index_str);
		flag = atoi(flag_str);

		// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
		//chr_index += (CODEBASE_START_BASE - SAM_START_BASE);
		_chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(flag & 0x10)
		{
			strand_char = 'R';
		}

		// Sanity check. Is this fragment mapped?
		if(flag & 0x04)
		{
			// The read is not mapping.
			chrom[0] = 0;
		}
		else
		{
			chr_index = _chr_index;
			sequenced_length = strlen(fragment);
		}
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_ELAND_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char cur_fragment[100];
	char quality_str[100];
	int _chr_index;
	char _strand_char;

	if(sscanf(cur_line, "%s %s %s %*d %*d %*d %s %d %c", read_id, cur_fragment, quality_str, chrom, &_chr_index, &_strand_char) == 6)
	{
		chr_index = _chr_index;
		sequenced_length = strlen(cur_fragment);
		sprintf(mapping_quality_str, "%dM", sequenced_length);

		strand_char = 'F';
		if(_strand_char == '-')
		{
			strand_char = 'R';
		}
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_bowtie_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	char nucs[1000];
	char strand_sign_str[100];
	char chr_start_index_str[100];
	if(sscanf(cur_line, "%[^'\t'] %[^'\t'] %[^'\t'] %[^'\t'] %[^'\t']", read_id, strand_sign_str, chrom, chr_start_index_str, nucs) == 5)
    {
		strand_sign = strand_sign_str[0];
		chr_start_index = atoi(chr_start_index_str);
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		chr_start_index = translate_coord(chr_start_index, BOWTIE_COORDS::start_base, CODEBASE_COORDS::start_base);

		sprintf(mapping_quality_str, "%dM", (int)strlen(nucs));

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
		sequenced_length = strlen(nucs);
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED4_read_line(char* cur_line, 
	char* read_id,
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;
	//if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
	if(sscanf(cur_line, "%s %d %d %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		sprintf(mapping_quality_str, "%dM", chr_end_index-chr_start_index);
		sequenced_length = chr_end_index-chr_start_index;

		chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
	}
	else
	{
		chrom[0] = 0;
	}
}

void preprocess_BED5_read_line(char* cur_line, 
	char* read_id, 
	char* chrom, 
	int& chr_index, int& sequenced_length, 
	char& strand_char, 
	char* mapping_quality_str)
{
	char strand_sign;
	int chr_start_index;
	int chr_end_index;
	if(sscanf(cur_line, "%s %d %d %*s %*s %c", chrom, &chr_start_index, &chr_end_index, &strand_sign) == 4)
    {
		// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
		sprintf(mapping_quality_str, "%dM", chr_end_index-chr_start_index);
		sequenced_length = chr_end_index-chr_start_index;

		chr_start_index = translate_coord(chr_start_index, BED_COORDS::start_base, CODEBASE_COORDS::start_base);

		// Check the flag and determine the strand.
		strand_char = 'F';
		if(strand_sign == '-')
		{
			strand_char = 'R';
		}

		chr_index = chr_start_index;
	}
	else
	{
		chrom[0] = 0;
	}
}

void delete_fragments(vector<t_mapped_fragment*>* fragment_list)
{
	for(int i_frag = 0; i_frag < (int)fragment_list->size(); i_frag++)
	{
		delete(fragment_list->at(i_frag));
	}

	fragment_list->clear();
	delete(fragment_list);
}

void delete_fragments(t_mapped_fragment** fragment_list)
{
	int i_frag = 0;
        while(fragment_list[i_frag] != NULL)
        {
                free(fragment_list[i_frag]);
		i_frag++;
        }

        delete[] fragment_list;
}

void load_fragments(char* mapped_reads_fp, 
	vector<t_mapped_fragment*>* fore_strand_frags, vector<t_mapped_fragment*>* rev_strand_frags, 
	int max_n_pcr_amplified_reads)
{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	fprintf(stderr, "Loading fragments from %s.\n", mapped_reads_fp);
}

	vector<t_mapped_read*>* fore_reads = new vector<t_mapped_read*>();
	vector<t_mapped_read*>* rev_reads = new vector<t_mapped_read*>();

	// Add the fragments: If the fragmens do not exist, the algorithm returns.
	//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
	load_reads(mapped_reads_fp, fore_reads, rev_reads, max_n_pcr_amplified_reads);
	get_mapped_fragments_per_mapped_reads(fore_reads, fore_strand_frags);
	get_mapped_fragments_per_mapped_reads(rev_reads, rev_strand_frags);

	// Delete the memory for the mapped reads, that is not needed any more.
	for(int i_f_r = 0; i_f_r < (int)fore_reads->size(); i_f_r++)
	{
		delete_mapped_read(fore_reads->at(i_f_r));
	} // i_f_r loop.
	delete fore_reads;

	for(int i_r_r = 0; i_r_r < (int)rev_reads->size(); i_r_r++)
	{
		delete_mapped_read(rev_reads->at(i_r_r));
	} // i_f_r loop.
	delete rev_reads;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	fprintf(stderr, "Loaded %ld forward and %ld reverse fragments.\n", fore_strand_frags->size(), rev_strand_frags->size());
}
}

//void load_reads_per_dir(char* mapped_reads_dir, vector<char*>* chr_ids, 
//	vector<vector<t_mapped_read*>*>* fore_strand_reads_per_chr, vector<vector<t_mapped_read*>*>* rev_strand_reads_per_chr, 
//	int max_n_pcr_amplified_reads)
//{
//	char cur_dir_chr_ids_fp[1000];
//	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
//	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);
//
//	if(cur_dir_chr_ids == NULL)
//	{
//		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
//		getc(stdin);
//		return;
//	}
//
//	// Load the fragment in this chromosome from all the directories.
//	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
//	{
//		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//
//		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));
//
//		char cur_dir_cur_chr_reads_fp[1000];
//		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));
//
//		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
//		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();
//
//		// Add the fragments: If the fragmens do not exist, the algorithm returns.
//		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
//		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);
//
//		// Add the loaded fragments.
//		fore_strand_reads_per_chr->push_back(cur_chr_fore_reads);
//		rev_strand_reads_per_chr->push_back(cur_chr_rev_reads);
//	} // i_chr loop.
//}

//void load_fragments_per_dir(char* mapped_reads_dir, 
//	vector<char*>* chr_ids, 
//	vector<vector<t_mapped_fragment*>*>* fore_strand_fragments_per_chr, 
//	vector<vector<t_mapped_fragment*>*>* rev_strand_fragments_per_chr, 
//	int max_n_pcr_amplified_reads)
//{
//	char cur_dir_chr_ids_fp[1000];
//	sprintf(cur_dir_chr_ids_fp, "%s/chr_ids.txt", mapped_reads_dir);
//	vector<char*>* cur_dir_chr_ids = buffer_file(cur_dir_chr_ids_fp);
//
//	if(cur_dir_chr_ids == NULL)
//	{
//		fprintf(stderr, "Could not load the chromosome id's from the directory %s\n", mapped_reads_dir);
//		getc(stdin);
//		return;
//	}
//
//	// Load the fragment in this chromosome from all the directories.
//	for(int i_chr = 0; i_chr < (int)cur_dir_chr_ids->size(); i_chr++)
//	{
//		//fore_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//		//rev_strand_fragments_per_chr->push_back(new vector<t_mapped_fragment*>());
//
//		chr_ids->push_back(t_string::copy_me_str(cur_dir_chr_ids->at(i_chr)));
//
//		char cur_dir_cur_chr_reads_fp[1000];
//		sprintf(cur_dir_cur_chr_reads_fp, "%s/%s_mapped_reads.txt", mapped_reads_dir, cur_dir_chr_ids->at(i_chr));
//
//		vector<t_mapped_read*>* cur_chr_fore_reads = new vector<t_mapped_read*>();
//		vector<t_mapped_read*>* cur_chr_rev_reads = new vector<t_mapped_read*>();
//
//		// Add the fragments: If the fragmens do not exist, the algorithm returns.
//		//load_fragments(cur_dir_cur_chr_reads_fp, fore_strand_fragments_per_chr->back(), rev_strand_fragments_per_chr->back(), max_n_pcr_amplified_reads, load_reads_only);
//		load_reads(cur_dir_cur_chr_reads_fp, cur_chr_fore_reads, cur_chr_rev_reads, max_n_pcr_amplified_reads);
//
//		vector<t_mapped_fragment*>* cur_chr_fore_frags = new vector<t_mapped_fragment*>(); 
//		get_mapped_fragments_per_mapped_reads(cur_chr_fore_reads, cur_chr_fore_frags);
//
//		vector<t_mapped_fragment*>* cur_chr_rev_frags = new vector<t_mapped_fragment*>(); 
//		get_mapped_fragments_per_mapped_reads(cur_chr_rev_reads, cur_chr_rev_frags);
//
//		// Add the loaded fragments.
//		fore_strand_fragments_per_chr->push_back(cur_chr_fore_frags);
//		rev_strand_fragments_per_chr->push_back(cur_chr_rev_frags);
//
//		// Delete the memory for the mapped reads, that is not needed any more.
//		delete_mapped_reads(cur_chr_fore_reads);
//		delete_mapped_reads(cur_chr_rev_reads);
//	} // i_chr loop.
//}

bool sort_mapped_fragments_per_3p(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	int frag1_3p = fragment_3p_accessor(frag1);
	int frag2_3p = fragment_3p_accessor(frag2);

	return(frag1_3p < frag2_3p);
}

bool sort_mapped_fragments(t_mapped_fragment* frag1, t_mapped_fragment* frag2)
{
	if(frag1->base_index < frag2->base_index)
	{
		return(frag1->base_index < frag2->base_index);
	}
	else if(frag1->base_index == frag2->base_index)
	{
		if((frag1->base_index + frag1->sequenced_fragment_length) < (frag2->base_index + frag2->sequenced_fragment_length))
		{
			return(true);
		}
		else
		{
			return(false);
		}
	}
	else
	{
		return(false);
	}
}

int fragment_3p_accessor(void* obj_ptr)
{
	t_mapped_fragment* frag_obj_ptr = (t_mapped_fragment*)obj_ptr;

	if(frag_obj_ptr->strand_char == 'F')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else if(frag_obj_ptr->strand_char == 'R')
	{
		return(frag_obj_ptr->base_index+frag_obj_ptr->sequenced_fragment_length-1);
	}
	else
	{
		fprintf(stderr, "The strand char for fragment object is %c @ %s(%d)\n", frag_obj_ptr->strand_char, __FILE__, __LINE__);
		exit(0);
		return(-1);
	}
}

vector<t_mapped_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_mapped_fragment*>* fore_frag_list, vector<t_mapped_fragment*>* rev_frag_list, int enrichment_mapped_fragment_length)
{
	vector<t_mapped_fragment*>* combined_frags = new vector<t_mapped_fragment*>();

	if(fore_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)fore_frag_list->size(); i_frag++)
		{
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));
			new_fragment_node->base_index = fore_frag_list->at(i_frag)->base_index;

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > fore_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				new_fragment_node->sequenced_fragment_length = fore_frag_list->at(i_frag)->sequenced_fragment_length;
			}


			new_fragment_node->strand_char = fore_frag_list->at(i_frag)->strand_char;

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// For the reverse fragment list, we need to also set the base index since it may be move further to left.
	if(rev_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < (int)rev_frag_list->size(); i_frag++)
		{
			// If this enrichment fragment length is set to 0, this means that there is no fragment extension.
			t_mapped_fragment* new_fragment_node = (t_mapped_fragment*)malloc(sizeof(t_mapped_fragment));

			// If the enrichment length is greater than the actual sequence tag length, update the sequenced length.
			if(enrichment_mapped_fragment_length > rev_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				// The enrichment fragment length is larger than the sequenced length, translate the fragment base.
				new_fragment_node->base_index = (rev_frag_list->at(i_frag)->base_index + rev_frag_list->at(i_frag)->sequenced_fragment_length - enrichment_mapped_fragment_length);
				new_fragment_node->sequenced_fragment_length = enrichment_mapped_fragment_length;
			}
			else
			{
				// The enrichment fragment length is smaller than the sequenced fragment length, do not do any translation of the fragment base.
				new_fragment_node->base_index = rev_frag_list->at(i_frag)->base_index;
				new_fragment_node->sequenced_fragment_length = rev_frag_list->at(i_frag)->sequenced_fragment_length;
			}

			new_fragment_node->strand_char = 'F'; // Change the place of these to forward strand.

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// Sort one last time.
	sort(combined_frags->begin(), combined_frags->end(), sort_mapped_fragments);

	return(combined_frags);
}

int* get_n_reads_per_window(int n_wins, vector<t_mapped_fragment*>* frag_list)
{
	int* n_reads_per_window = new int[n_wins + 2];
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		n_reads_per_window[i_win] = 0;
	} // i_win loop

	for(int i_frag = 0; i_frag < (int)frag_list->size(); i_frag++)
	{
		int cur_i_win = frag_list->at(i_frag)->base_index / MEG_BASE;

		n_reads_per_window[cur_i_win]++;
	} // i_nuc loop.

	return(n_reads_per_window);
}

enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
	int i = 0;

	is_read_spliced = false;

	int state = VAL;
	while(mapping_map_str[i] != 0)
	{		
		if(state == VAL)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
			// MIDNSHPX=
			if(mapping_map_str[i] == 'M' ||
				mapping_map_str[i] == 'I' ||
				mapping_map_str[i] == 'D' ||
				mapping_map_str[i] == 'N' ||
				mapping_map_str[i] == 'S' ||		
				mapping_map_str[i] == 'H' ||
				mapping_map_str[i] == 'P' ||
				mapping_map_str[i] == 'X' ||
				mapping_map_str[i] == '=')
			{
				state = TYPE;

				if(mapping_map_str[i] != 'M')
				{
					is_read_spliced = true;
				}
			}
			else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				// State is still VAL.
			}
			else
			{
				return(false);
			}
		}
		else if(state == TYPE)
		{
			//fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
			// A number is expected.
			if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
			{
				state = VAL;
			}
			else
			{
				return(false);
			}
		}

		// Move to next character.
		i++;
	}

	return(true);
}

void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										//t_string* cur_entry_length_str,
										int& l_cur_entry,
										char& entry_type_char)
{	
	// Clean the length string.
	//cur_entry_length_str->empty();
	l_cur_entry = 0;

	// Get the next entry in the cigar string.
	while(mapping_map_str[i_mapp_map] != 0)
	{
		if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
		{
			break;
		}
		//cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
		l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
		i_mapp_map++;
	}

	is_matching = false;
	if(mapping_map_str[i_mapp_map] == 'M')
	{
		//fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
		is_matching = true;
	}
	else
	{
		//fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
	}	

	entry_type_char = mapping_map_str[i_mapp_map];

	// Move over the current entry identifier.
	i_mapp_map++;
}

/*
Load reads from a preprocesed read file.
*/
void load_reads(char* mapped_reads_fp, 
	vector<t_mapped_read*>* pruned_fore_strand_reads, vector<t_mapped_read*>* pruned_rev_strand_reads, 
	int max_n_pcr_amplified_reads)
{
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	printf("Loading mapped-reads from %s.\n", mapped_reads_fp);

	if(!check_file(mapped_reads_fp))
	{
		printf("Could not open mapped reads file %s @ %s(%d).\n", mapped_reads_fp, __FILE__, __LINE__);
		return;
	}

	vector<t_mapped_read*>* fore_strand_reads = new vector<t_mapped_read*>();
	vector<t_mapped_read*>* rev_strand_reads = new vector<t_mapped_read*>(); 

	//char cur_fragment[10000];
	char mapping_map_str[10000];
	char strand_char;
	int chr_index;

	// Buffer the whole file.
	t_file_buffer* mapped_read_file_buffer = load_file(mapped_reads_fp);

	// Read and validate the mapped reads in the file.
	//while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	while(1)
	{
		char* cur_line = getline_per_file_buffer(mapped_read_file_buffer);

		if(cur_line == NULL)
		{
			break;
		}

		if(sscanf(cur_line, "%s %c %d", mapping_map_str, &strand_char, &chr_index) != 3)
		{
			fprintf(stderr, "Could not parse fragment line: %s\n", cur_line);
		}

		int i_mapp_map = 0;
		t_string* cur_entry_length_str = new t_string();
		bool is_matching = false;
		char entry_type_char;

		//fprintf(stderr, "Processing cigar string: %s\n", quality_str);
		bool is_read_spliced = false;
		bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

		// Allocate and initialize the new mapped read.
		t_mapped_read* new_mapped_read = new t_mapped_read();
		new_mapped_read->mapping_str = t_string::copy_me_str(mapping_map_str);
		new_mapped_read->strand = strand_char;
		new_mapped_read->span = 0;
		int left_posn = chr_index;

		// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
		while(mapping_map_str_valid && 
			mapping_map_str[i_mapp_map] != 0)
		{
			int l_cur_entry;
			get_next_entry_per_mapp_map_string(mapping_map_str,
												i_mapp_map, 
												is_matching,
												l_cur_entry,
												entry_type_char);

			// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.

			// Update the base for the current entry.
			// Must check whether to update the mapping posn: Update only for D and M entries.
			/*if(entry_type_char == 'D' || 
				entry_type_char == 'M' ||
				entry_type_char == 'N' ||
				entry_type_char == 'H')*/
			if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
			{
				chr_index += l_cur_entry;
				new_mapped_read->span += l_cur_entry;
			}
		} // mapping map string processing loop.

		// Set the base_index.
		new_mapped_read->base_index = left_posn;	

		new_mapped_read->strand = new_mapped_read->strand;

		// Add to the lists.
		if(new_mapped_read->strand == 'F')
		{
			fore_strand_reads->push_back(new_mapped_read);
		}
		else if(new_mapped_read->strand == 'R')
		{
			rev_strand_reads->push_back(new_mapped_read);
		}

		delete(cur_entry_length_str);
		delete [] cur_line;
	} // current fragment data reading loop.

	// Free file buffer.
	unload_file(mapped_read_file_buffer);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
	printf("Loaded %ld fragments on forward strand.\n", fore_strand_reads->size());
    printf("Loaded %ld fragments on reverse strand.\n", rev_strand_reads->size());
}

	// Prune the reads.
	vector<t_mapped_read*>* all_reads = new vector<t_mapped_read*>();
	all_reads->insert(all_reads->end(), fore_strand_reads->begin(), fore_strand_reads->end());
	all_reads->insert(all_reads->end(), rev_strand_reads->begin(), rev_strand_reads->end());
	
	// Following gets rid of the pruned reads and returns only the pruned reads without re-allocating the reads/split-mapped-fragments.
	prune_reads(all_reads, max_n_pcr_amplified_reads, 
				pruned_fore_strand_reads, 
				pruned_rev_strand_reads);

	// Can clean the memory for all reads. Note that the pruned reads are cleaned up above, unpruned ones are put in the reads per strand lists.
	delete(fore_strand_reads);
	delete(rev_strand_reads);
	delete(all_reads);
}

void add_mapped_fragments_per_mapped_read(t_mapped_read* mapped_read, vector<t_mapped_fragment*>* mapped_fragments)
{
	int i_mapp_map = 0;
	t_string* cur_entry_length_str = new t_string();
	bool is_matching = false;
	char entry_type_char;
	char strand_char = mapped_read->strand;
	//int chr_index = (mapped_read->strand=='F')?(mapped_read->base_index):(mapped_read->base_index-mapped_read->span+1);
	int chr_index = (mapped_read->base_index);
	char* mapping_map_str = mapped_read->mapping_str;

	//fprintf(stderr, "Processing cigar string: %s (%d, %c)\n", mapping_map_str, chr_index, strand_char);
	bool is_read_spliced = false;
	bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

	// If loading of all the fragments requested, process all the entries in the mapping string, then add all of them to the list of fragments for this read.
	while(mapping_map_str_valid && 
		mapping_map_str[i_mapp_map] != 0)
	{
		int l_cur_entry = 0;
		get_next_entry_per_mapp_map_string(mapping_map_str,
											i_mapp_map, 
											is_matching,
											l_cur_entry,
											entry_type_char);

		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.		

		//int l_fragment = strlen(cur_fragment);
		if(is_matching)
		{
			//int l_fragment = get_l_fragment_per_cigar(quality_str);
			if(strand_char == 'F')
			{
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;
		
				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
			}
			else if(strand_char == 'R')
			{
				// Allocate and initialize a fragment and add it to the reverse strand fragment list.			
				t_mapped_fragment* new_fragment = new t_mapped_fragment();
				//new_fragment->base_index = chr_index + l_cur_entry - 1;
				new_fragment->base_index = chr_index;
				new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = l_cur_entry;

				mapped_fragments->push_back(new_fragment);

				//fprintf(stderr, "Adding: %c, %d (%d)\n", new_fragment->strand_char, new_fragment->base_index, new_fragment->sequenced_fragment_length);
				//getc(stdin);
				//rev_strand_frags->push_back(new_fragment);
			} // reverse strand check.
		} // maching check.

		// Update the base for the current entry.
		// Must check whether to update the mapping posn: Update only for D and M entries.
		//if(entry_type_char == 'D' || 
		//	entry_type_char == 'M' ||
		//	entry_type_char == 'N' ||
		//	entry_type_char == 'H')
		if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
		{
			chr_index += l_cur_entry;
		}
	} // mapping map string processing loop.

	delete cur_entry_length_str;
}

// Check if the 
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'H' ||
	//if(entry_char == 'P' ||
	//if(entry_char == 'X' ||
	if(entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}

bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'D' || 
	//	entry_char == 'M' ||
	//	entry_char == 'N' ||
	//	entry_char == 'H')
	if(entry_char == 'D' || 
		entry_char == 'M' ||
		entry_char == 'N')
	{
		return(true);
	}

	return(false);
}

void get_mapped_fragments_per_mapped_reads(vector<t_mapped_read*>* mapped_reads, vector<t_mapped_fragment*>* mapped_fragments)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		add_mapped_fragments_per_mapped_read(mapped_reads->at(i_r), mapped_fragments);
	} // i_r loop.
}

void delete_mapped_reads(vector<t_mapped_read*>* mapped_reads)
{
	for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
	{
		delete_mapped_read(mapped_reads->at(i_r));
	} // i_r loop.

	delete(mapped_reads);
}

void delete_mapped_read(t_mapped_read* mapped_read)
{
	delete [] mapped_read->mapping_str;
	delete(mapped_read);
}

/*
Prune reads:
Note that the trick here is to deal with the reads that have exact same pattern of mapping. We do not care about the
strand, this should be taken care of before the function is called.

Note that the pruning must be done at the read level, not at the fragment level.
*/
void prune_reads(vector<t_mapped_read*>* mapped_reads, int n_max_reps_per_posn, 
	vector<t_mapped_read*>* pruned_forward_reads, 
	vector<t_mapped_read*>* pruned_reverse_reads)
{
	// If the pruning is not requested, return all the reads.
	if(n_max_reps_per_posn == 0)
	{
		fprintf(stderr, "Skipping pruning.\n");
		for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		} // i_r loop.

		return;
	}

	// Sort the mapped reads with respect to their 5' posn.
	sort(mapped_reads->begin(), mapped_reads->end(), sort_mapped_reads_per_5p);

    // First get rid of the extra fragments on forward strand.
	int* rep_cnts = new int[mapped_reads->size() + 2];
	memset(rep_cnts, 0, sizeof(int) * (mapped_reads->size() + 1));

	int prev_read_left_posn = 0;
    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		int cur_read_left_posn = (mapped_reads->at(i_r)->base_index);
		
		if(i_r > 0 &&
			prev_read_left_posn == cur_read_left_posn)
		{
			rep_cnts[i_r] = rep_cnts[i_r-1] + 1;
		}
		else // This is a new fragment set its copy number to 1.
		{
			rep_cnts[i_r] = 1;
		}

		// Update prev. left posn.
		prev_read_left_posn = cur_read_left_posn;
    } // i_frag loop.

    for(int i_r = 0; i_r < (int)mapped_reads->size(); i_r++)
    {
		if(rep_cnts[i_r] <= n_max_reps_per_posn)
		{
			if(mapped_reads->at(i_r)->strand == 'F')
			{
				pruned_forward_reads->push_back(mapped_reads->at(i_r));
			}
			else
			{
				pruned_reverse_reads->push_back(mapped_reads->at(i_r));
			}
		}
		else // This is a new fragment set its copy number to 1.
		{
			// Delete the pruned reads, otherwise they will be lost.
			delete_mapped_read(mapped_reads->at(i_r));

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
			fprintf(stderr, "Pruning %d repetition read @ %d, %c.\n", rep_cnts[i_r], mapped_reads->at(i_r)->base_index, mapped_reads->at(i_r)->strand);
			getc(stdin);
}
		}
    } // i_r loop.

	delete [] rep_cnts;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
	fprintf(stderr, "Pruned to %ld forward, %ld reverse strand fragments.\n", pruned_forward_reads->size(), pruned_reverse_reads->size());
}

bool sort_mapped_reads_per_5p(t_mapped_read* read1, t_mapped_read* read2)
{
	int frag1_5p = read_5p_accessor(read1);
	int frag2_5p = read_5p_accessor(read2);

	return(frag1_5p < frag2_5p);
}

int read_5p_accessor(void* obj_ptr)
{
	t_mapped_read* frag_obj_ptr = (t_mapped_read*)obj_ptr;

	return(frag_obj_ptr->base_index);
}


