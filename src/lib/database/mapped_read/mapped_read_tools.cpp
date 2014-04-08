#include <stdio.h>
#include <stdlib.h>
#include "mapped_read_tools.h"
#include <vector>
#include <ctype.h>
#include "../../../lib/chromosome/chromosome.h"
#include "../../../lib/database/signal_track/signal_track_tools.h"
#include "../../../lib/chromosome/indexing_per_format.h"
#include "../../../lib/utils/file/utils.h"
#include "fragment.h"
#include <string.h>


using namespace std; 

/*
Mapped read file interface:
===========================
Contains the functions for reformatting the mapped read files in different formats to text based format that can be loaded using the interface defined in fragment.h and also for writing the binary based fragment file. 
This interface converts all the different formats of mapped tag data to a single format that is used by PeakSeq. There are some conventions while reformatting additional file formats:
- All the indices are converted to 1 based indices. The fragment.cpp file assumes indices are 1 based while validating the reads. Also chromosome class expects 1 based indices when accessing the nucleotides.
- The strand that gets dumped to the text based output file for reverse strands have to be the strand on the reverse strand. In other words, the tags in the
text based file are the tags themselves. For example SAM files contain the reverse-complement of the tag for the tags that map to negative strand. Those 
should be reverse complemented back to the tag itself while writing the text based op file.
- The indexing conversions must be handled correctly, some file formats store indices as 0-based and some as 1-based.

The delta profile site list corresponds to 
*/

FILE* get_fragment_file_pointer_by_chr_file_name(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_fn)
{
        for(int i = 0; i < chr_ids->size(); i++)
        {
			// Note that the id is not the file name. To get the file name, must add .fa to the id.
			char id_fn[1000];
			sprintf(id_fn, "%s.fa", chr_ids->at(i));
			if(strcmp(id_fn, chr_fn) == 0)
			{
					return(fragment_f_ptrs->at(i));
			}
        } // i loop

        return(NULL);
}

FILE* get_fragment_file_pointer_by_chr_id(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_id)
{
        for(int i = 0; i < chr_ids->size(); i++)
        {
		// Note that the id is not the file name. To get the file name, must add .fa to the id.
                if(strcmp(chr_ids->at(i), chr_id) == 0)
                {
                        return(fragment_f_ptrs->at(i));
                }
        } // i loop

        return(NULL);
}

void parse_ELAND_mapped_reads_file(vector<char*>* chr_ids, char* parsed_reads_op_dir, char* eland_fp)
{
	fprintf(stderr, "Parsing ELAND formatted mapped reads file %s\n", eland_fp);

        // Divide ELAND output with respect to chromosomes.
		FILE* f_eland = NULL;
		if(strcmp(eland_fp, "stdin") == 0)
		{
			f_eland = stdin;
		}
		else
		{
			f_eland = open_f(eland_fp, "r");
		}

        char cur_line[2000];
        int n_frags = 0;
        int n_total_frags = 0;

        vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
        for(int i = 0; i < chr_ids->size(); i++)
        {
                char new_fn[10000];

				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i));
	             frag_f_ptrs->push_back(open_f(new_fn, "w"));
        }

        while(1)
        {
                if(fgets(cur_line, 1000, f_eland) != NULL)
                {
                        // >FC20B5RA_20080331:1:1:121:603  TTTTCAAAAGTAGTGGTTCTTTCCTTAT    U0      1       0       0       chr13.fa 46121449       F       ..
                        char cur_fragment[100];
                        char quality_str[100];
                        char chr_fn[300];
                        int chr_index;
                        char strand_char;
                        //fprintf(stderr, "Read: %s\n", cur_line);

                        n_total_frags++;

                        if(sscanf(cur_line, "%*s %s %s %*d %*d %*d %s %d %c", cur_fragment, quality_str, chr_fn, &chr_index, &strand_char) == 5)
                        {
								chr_index += (PEAKSEQ_BASE - ELAND_BASE); // change the index.

                                //fprintf(stderr, "%s [%c:%d @ %s]     \r", cur_fragment, strand_char, chr_index, quality_str);

                                // Extend and map the position on the forward chromosome coordinates.
                                //if(n_frags % 1000 == 0)
                                //{
                                        //fprintf(stderr, "%d. fragment.          \r", n_frags);
                                //}

                                // Open the chromosome file output for current entry.
                                FILE* cur_frag_file = get_fragment_file_pointer_by_chr_id(frag_f_ptrs, chr_ids, chr_fn);
                                if(cur_frag_file == NULL)
                                {
                                        printf("Could not resolve file pointer for fragment with file name %s\n", chr_fn);
                                }
								else
								{
									fprintf(cur_frag_file, "%s %s %c %d\n", cur_fragment, quality_str, strand_char, chr_index);
									n_frags++;
								}
                        }

                        if(n_total_frags % 100000 == 0)
                        {
                                fprintf(stderr, "%d (%d)         \r", n_total_frags, n_frags);
                        }
                }
                else
                {
                        break;
                }
        } // File parsing loop.i

        fclose(f_eland);

        for(int i = 0; i < frag_f_ptrs->size(); i++)
        {
                fclose(frag_f_ptrs->at(i));
        }

        frag_f_ptrs->clear();
        delete(frag_f_ptrs);
}

void parse_tagAlign_formatted_mapped_reads_file(vector<char*>* chr_ids, char* parsed_reads_op_dir, char* tagAlign_fp)
{
		fprintf(stderr, "Parsing tagAlign formatted mapped reads file %s\n", tagAlign_fp);

       // Divide SAM output with respect to chromosomes.
        FILE* f_tagAlign = NULL;
		if(strcmp(tagAlign_fp, "stdin") == 0)
		{
			f_tagAlign = stdin;
		}
		else
		{
			f_tagAlign = open_f(tagAlign_fp, "r");
			//if(f_tagAlign == NULL)
			//{
			//		printf("Could not open %s\n", tagAlign_fp);
			//		exit(0);
			//}
		}

        char cur_line[2000];
        int n_frags = 0;
        int n_total_frags = 0;

        vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
        for(int i = 0; i < chr_ids->size(); i++)
        {
                char new_fn[10000];

				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i));
	            frag_f_ptrs->push_back(open_f(new_fn, "w"));
        }

        while(1)
        {
                if(fgets(cur_line, 1000, f_tagAlign) != NULL)
                {
						//chrX 8823384 8823409 AGAAGGAAAATGATGTGAAGACATA 1000 +
                        char quality_str[100];
                        char chr_id[300];
                        int chr_start_index;
						int chr_end_index;
                        char strand_sign;
                        //fprintf(stderr, "Read: %s\n", cur_line);

                        n_total_frags++;

                        //if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %*s", &flag, chr_id, &chr_index, cigar_str, cur_fragment) == 5)
						if(sscanf(cur_line, "%s %d %d %*s %*d %c", chr_id, &chr_start_index, &chr_end_index, &strand_sign) == 4)
                        {
							// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
							chr_start_index += (PEAKSEQ_BASE - tagAlign_BASE);
							chr_end_index += (PEAKSEQ_BASE - tagAlign_BASE);

							strcpy(quality_str, "XX");

							// Check the flag and determine the strand.
							char strand_char = 'F';
							if(strand_sign == '-')
							{
								strand_char = 'R';
							}

							//fprintf(stderr, "%s [%c:%d @ %s]     \r", cur_fragment, strand_char, chr_index, quality_str);

							// Extend and map the position on the forward chromosome coordinates.
							//if(n_frags % 1000 == 0)
							//{
									//fprintf(stderr, "%d. fragment.          \r", n_frags);
							//}

							// Open the chromosome file output for current entry.
							FILE* cur_frag_file = get_fragment_file_pointer_by_chr_id(frag_f_ptrs, chr_ids, chr_id);

							// Check if the current entry has a file pointer associated with it, if not, skip it.
							
							if(cur_frag_file == NULL)
							{
									printf("Could not resolve file pointer for chromosome id %s\n", chr_id);
							}
							else
							{								
								// The fragment information is skipped in the tagAlign files.
								char* cur_fragment = new char[chr_end_index - chr_start_index + 5];
								int n_frag_nucs = chr_end_index - chr_start_index + 1;
								for(int i_nuc = 0; i_nuc < n_frag_nucs; i_nuc++)
								{
									cur_fragment[i_nuc] = 'N';
								} // i_nuc loop.
								cur_fragment[n_frag_nucs] = 0;

								fprintf(cur_frag_file, "%s %s %c %d\n", cur_fragment, quality_str, strand_char, chr_start_index);
								n_frags++;
							}
                        }

                        if(n_total_frags % 100000 == 0)
                        {
                                fprintf(stderr, "%d (%d)         \r", n_total_frags, n_frags);
                        }
                }
                else
                {
                        break;
                }
        } // File parsing loop.i

        fclose(f_tagAlign);

        for(int i = 0; i < frag_f_ptrs->size(); i++)
        {
                fclose(frag_f_ptrs->at(i));
        }

        frag_f_ptrs->clear();
        delete(frag_f_ptrs);
}

void parse_bowtie_formatted_mapped_reads_file(vector<char*>* chr_ids, char* parsed_reads_op_dir, char* bowtie_fp)
{
		fprintf(stderr, "Parsing bowtie formatted mapped reads file %s\n", bowtie_fp);

       // Divide SAM output with respect to chromosomes.
        FILE* f_bowtie = NULL;
		if(strcmp(bowtie_fp, "stdin") == 0)
		{
			f_bowtie = stdin;
		}
		else
		{
			f_bowtie = open_f(bowtie_fp, "r");
			//if(f_tagAlign == NULL)
			//{
			//		printf("Could not open %s\n", tagAlign_fp);
			//		exit(0);
			//}
		}

        char cur_line[2000];
        int n_frags = 0;
        int n_total_frags = 0;

        vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
        for(int i = 0; i < chr_ids->size(); i++)
        {
                char new_fn[10000];

				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i));
	            frag_f_ptrs->push_back(open_f(new_fn, "w"));
        }

        while(1)
        {
                if(fgets(cur_line, 1000, f_bowtie) != NULL)
                {
						//2:SOLEXA-1GA-2:2:97:1089:2007#0 +       chr1_maternal   242     ACCCTAAACCCTAAACCCTAACCCTAACCCTAACCC    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    0
                        char nucs[1000];
						char quality_str[10];
                        char chr_id[300];
                        int chr_start_index;
						int chr_end_index;
                        char strand_sign;
                        //fprintf(stderr, "Read: %s\n", cur_line);

                        n_total_frags++;

                        //if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %*s", &flag, chr_id, &chr_index, cigar_str, cur_fragment) == 5)
						//if(sscanf(cur_line, "%s %d %d %*s %*d %c", chr_id, &chr_start_index, &chr_end_index, &strand_sign) == 4)
						if(sscanf(cur_line, "%*s %c %s %d %s", &strand_sign, chr_id, &chr_start_index, nucs) == 4)
                        {
							// Note that the indices in tagAlign file are 0 based, these must be translated to 1 based indices.
							chr_start_index += (PEAKSEQ_BASE - bowtie_BASE);
							//chr_end_index += (PEAKSEQ_BASE - bowtie_BASE);

							strcpy(quality_str, "XX");

							// Check the flag and determine the strand.
							char strand_char = 'F';
							if(strand_sign == '-')
							{
								strand_char = 'R';
							}

							// Open the chromosome file output for current entry.
							FILE* cur_frag_file = get_fragment_file_pointer_by_chr_id(frag_f_ptrs, chr_ids, chr_id);

							// Check if the current entry has a file pointer associated with it, if not, skip it.
							
							if(cur_frag_file == NULL)
							{
									printf("Could not resolve file pointer for chromosome id %s\n", chr_id);
							}
							else
							{								
								fprintf(cur_frag_file, "%s %s %c %d\n", nucs, quality_str, strand_char, chr_start_index);
								n_frags++;
							}
                        }

                        if(n_total_frags % 100000 == 0)
                        {
                                fprintf(stderr, "%d (%d)         \r", n_total_frags, n_frags);
                        }
                }
                else
                {
                        break;
                }
        } // File parsing loop.i

        fclose(f_bowtie);

        for(int i = 0; i < frag_f_ptrs->size(); i++)
        {
                fclose(frag_f_ptrs->at(i));
        }

        frag_f_ptrs->clear();
        delete(frag_f_ptrs);
}

void parse_SAM_formatted_mapped_reads_file(vector<char*>* chr_ids, char* parsed_reads_op_dir, char* sam_fp)
{
		fprintf(stderr, "Parsing SAM formatted mapped reads file %s\n", sam_fp);

        // Divide SAM output with respect to chromosomes.
        FILE* f_sam = NULL;
		if(strcmp(sam_fp, "stdin") == 0)
		{
			f_sam = stdin;
		}
		else
		{
			f_sam = open_f(sam_fp, "r");
			//if(f_sam == NULL)
			//{
			//		printf("Could not open %s\n", sam_fp);
			//		exit(0);
			//}
		}

        char cur_line[2000];
        int n_frags = 0;
        int n_total_frags = 0;

        vector<FILE*>* frag_f_ptrs = new vector<FILE*>();
        for(int i = 0; i < chr_ids->size(); i++)
        {
                char new_fn[10000];

				sprintf(new_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i));
	             frag_f_ptrs->push_back(open_f(new_fn, "w"));
        }

        while(1)
        {
                if(fgets(cur_line, 1000, f_sam) != NULL)
                {
						// SOLEXA-1GA-2_0001:6:49:1106:529#0       16      chr1    10149   25      36M     *       0       0       CCTACCCCTAACCCTAACCCTAACCCTAACCTAACC    Z^[^IRRRH`S`ZWZba``ZJZ[`bbb`]bbaaaJT       NM:i:1  X1:i:1  MD:Z:31C4
                        char cur_fragment[100];
                        char quality_str[100];
                        char chr_id[300];
                        int chr_index;
                        char strand_char;
						char cigar_str[100];
						int flag;
                        //fprintf(stderr, "Read: %s\n", cur_line);

                        n_total_frags++;

                        if(sscanf(cur_line, "%*s %d %s %d %*d %s %*s %*d %*d %s %*s", &flag, chr_id, &chr_index, cigar_str, cur_fragment) == 5)
                        {
								// Translate the 0 based index in SAM file to ELAND's 1 based indexing.
								chr_index += (PEAKSEQ_BASE - SAM_BASE);;

								strcpy(quality_str, cigar_str);

								// Check the flag and determine the strand.
								strand_char = 'F';
								if(flag & 0x10)
								{
									strand_char = 'R';
								}

								// Sanity check. Is this fragment mapped?
								if(flag & 0x04)
								{
									// Quietly skip the tags that were not mapped.
									//printf("This fragment is unmapped what should i do:\n%s\n", cur_line);
									//getc(stdin);
								}
								else
								{
									//fprintf(stderr, "%s [%c:%d @ %s]     \r", cur_fragment, strand_char, chr_index, quality_str);

									// Extend and map the position on the forward chromosome coordinates.
									//if(n_frags % 1000 == 0)
									//{
											//fprintf(stderr, "%d. fragment.          \r", n_frags);
									//}

									// Open the chromosome file output for current entry.
									FILE* cur_frag_file = get_fragment_file_pointer_by_chr_id(frag_f_ptrs, chr_ids, chr_id);

									// Check if the current entry has a file pointer associated with it, if not, skip it.
									if(cur_frag_file == NULL)
									{
											printf("Could not resolve file pointer for chromosome id %s\n", chr_id);
									}
									else
									{
										// Reverse complement the current fragment.
										if(strand_char == 'R')
										{
											char* rev_comp_frag = new char[strlen(cur_fragment) + 2];
											memset(rev_comp_frag, 0, sizeof(char) * (strlen(cur_fragment) + 1));
	
											char rev_nucs[] = "NACGTN";
											int n_frag_nucs =  strlen(cur_fragment);
											for(int i_nuc = 0; i_nuc < n_frag_nucs; i_nuc++)
											{
												if(toupper(cur_fragment[i_nuc]) == 'A')
												{
													rev_comp_frag[n_frag_nucs - 1 - i_nuc] = 'T';
												}
												else if(toupper(cur_fragment[i_nuc]) == 'C')
                                                {
                                                        rev_comp_frag[n_frag_nucs - 1 - i_nuc] = 'G';
                                                }
                                                else if(toupper(cur_fragment[i_nuc]) == 'G')
                                                {
                                                        rev_comp_frag[n_frag_nucs - 1 - i_nuc] = 'C';
                                                }
                                                else if(toupper(cur_fragment[i_nuc]) == 'T')
                                                {
                                                        rev_comp_frag[n_frag_nucs - 1 - i_nuc] = 'A';
                                                }
												else
												{
													rev_comp_frag[n_frag_nucs - 1 - i_nuc] = 'N';
												}
											} // i_nuc loop.

											//printf("%s\n%s\n", cur_fragment, rev_comp_frag);
											//getc(stdin);
											strcpy(cur_fragment, rev_comp_frag);
											delete [] rev_comp_frag;
										}

										fprintf(cur_frag_file, "%s %s %c %d\n", cur_fragment, quality_str, strand_char, chr_index);
										n_frags++;
									}
								} // Mapped check for fragment.
                        }

                        if(n_total_frags % 100000 == 0)
                        {
                                fprintf(stderr, "%d (%d)         \r", n_total_frags, n_frags);
                        }
                }
                else
                {
                        break;
                }
        } // File parsing loop.i

        fclose(f_sam);

        for(int i = 0; i < frag_f_ptrs->size(); i++)
        {
                fclose(frag_f_ptrs->at(i));
        }

        frag_f_ptrs->clear();
        delete(frag_f_ptrs);
}

/*
Uses the smae directory to dump the binary files.
*/
void validate_dump_binary_files(vector<char*>* chr_ids, char* parsed_reads_op_dir, bool validate)
{
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
        char mapped_reads_fn[10000];
		char mapped_reads_bin_fn[10000];

        sprintf(mapped_reads_fn, "%s/%s_mapped_reads.txt", parsed_reads_op_dir, chr_ids->at(i_chr));
        sprintf(mapped_reads_bin_fn, "%s/%s_mapped_reads.bin", parsed_reads_op_dir, chr_ids->at(i_chr));

		char cur_chr_fp[10000];
		sprintf(cur_chr_fp, "%s.fa", chr_ids->at(i_chr));
		printf("Validating the reads in %s\n", cur_chr_fp);

		// Load the current chromosome, if validation is requested.
		t_chromosome* cur_chr = NULL;
		if(validate)
		{
			cur_chr = new t_chromosome(cur_chr_fp);
		}

		// Load the fragments from the text file for current chromosome.
        vector<t_fragment*>* fore_strand_frags = new vector<t_fragment*>();
        vector<t_fragment*>* rev_strand_frags = new vector<t_fragment*>();
		
		// Load the fragment lists.
        load_fragments(cur_chr, mapped_reads_fn, fore_strand_frags, rev_strand_frags);

		// Concatenate the lists of fragments.
		vector<t_fragment*>* all_frags = new vector<t_fragment*>();
		for(int i_frag = 0; i_frag < fore_strand_frags->size(); i_frag++)
		{
			all_frags->push_back(fore_strand_frags->at(i_frag));
		}

        for(int i_frag = 0; i_frag < rev_strand_frags->size(); i_frag++)
        {
                all_frags->push_back(rev_strand_frags->at(i_frag));			
        }

		// Dump the binary fragments.
		dump_fragments_binary(mapped_reads_bin_fn, all_frags);

		delete_fragments(fore_strand_frags);
		delete_fragments(rev_strand_frags);

		// Careful while deleting memory for all fragments.
		all_frags->clear();
		delete(all_frags);

		if(cur_chr != NULL)
		{
			delete(cur_chr);
		}
	} // i_chr loop.
}

double get_n_mapped_nucs(vector<t_fragment*>* fragments)
{
	double n_mapped_nucs = 0.0;

	for(int i_frag = 0; i_frag < fragments->size(); i_frag++)
	{
		n_mapped_nucs += fragments->at(i_frag)->sequenced_fragment_length;
	} // i_frag loop.

	return(n_mapped_nucs);
}

/*
Normalizes the profile with respect to number of mapped nucleotides in million, this is the only function that does the 
normalization since this is the only place where the number of mapped nucleotides can be computed. The bedGraphs do not 
include this information, thus signal_track_tools library will include, presumably, the normalized signal data.
*/
vector<t_profile_site*>* buffer_normalized_profile(vector<t_fragment*>* fragments, int min_i_nuc, int max_i_nuc, bool normalize)
{
	int height = 0;

	double n_total_mapped_nucs = get_n_mapped_nucs(fragments);

	double n_total_mapped_nucs_in_million = n_total_mapped_nucs / (1000 * 1000);

	// Generate the fragment ends.
	vector<t_frag_end*>* frag_ends = get_frag_ends(fragments);

	vector<t_profile_site*>* profile_sites = new vector<t_profile_site*>();

	int prev_pos = 0;
	int prev_height = 0;
	int cur_pos = 1;
    for(int i_end = 0; i_end < frag_ends->size(); i_end++)
    {
		/*
		A new region is added between prev_pos and this->frag_ends->at(i_end). If the requested region overlaps
		with the added region, dump the previous position. Note that the height in the added region is the hiehgt
		that is updated at the previous position.
		*/
        // Update the current height: Note that there can be multiple ends at this position.
        height += frag_ends->at(i_end)->side;
        cur_pos = frag_ends->at(i_end)->pos;

        // Now add all the other ends at this position.
        while((i_end + 1) < frag_ends->size() &&
                    cur_pos == frag_ends->at(i_end+1)->pos)
        {
                i_end++;
                height += frag_ends->at(i_end)->side;
        }

		if(prev_pos >= max_i_nuc)
		{
			// The profile passed the requested region completely, can return from the function, dump the previous position and height and return, this is necessary to account for the last position.
                        //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);

			t_profile_site* new_site = new t_profile_site();
			new_site->i_nuc = prev_pos;
			new_site->height = prev_height;
			profile_sites->push_back(new_site);

			break;
		}
		else if(cur_pos < min_i_nuc)
		{
			// Still did not reach the requested region.
		}
		else
		{
			// The newly added region overlaps with the requested region, dump the previous region and height, since 
			// the region corresponds to the height of the revious pos.
            //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);
			t_profile_site* new_site = new t_profile_site();
			new_site->i_nuc = prev_pos;
			new_site->height = prev_height;
			profile_sites->push_back(new_site);
		}
	
		// Update the previous position before moving to the next region.
		prev_pos = cur_pos;
		prev_height = height;
/*
        if(this->frag_ends->at(i_end)->pos >= min_i_nuc &&
                this->frag_ends->at(i_end)->pos < max_i_nuc)
        {
			fprintf(f_prof, "%d %d\n", this->frag_ends->at(i_end)->pos, height);
        }
*/
    } // i_frag

	// Normalize with respect to the total # of mapped nucleotides, if it is requested.
	if(normalize)
	{
		for(int i_site = 0; i_site < profile_sites->size(); i_site++)
		{
			profile_sites->at(i_site)->height /= n_total_mapped_nucs_in_million;
		} // i_site loop.
	}

	// Clean memory.
	for(int i_frag = 0; i_frag < frag_ends->size(); i_frag++)
	{
		delete(frag_ends->at(i_frag));
	} // i_frag loop.

	delete(frag_ends);

	return(profile_sites);
}

