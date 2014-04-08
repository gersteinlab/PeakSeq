#include <stdio.h>
#include <stdlib.h>
#include "fragment.h"
#include "../../../lib/chromosome/chromosome.h"
#include <algorithm>
#include <ctype.h>
#include "../../../lib/utils/file/utils.h"
#include <string.h>

void delete_fragments(vector<t_fragment*>* fragment_list)
{
	for(int i_frag = 0; i_frag < fragment_list->size(); i_frag++)
	{
		free(fragment_list->at(i_frag));
	}

	fragment_list->clear();
	delete(fragment_list);
}

void delete_fragments(t_fragment** fragment_list)
{
	int i_frag = 0;
        while(fragment_list[i_frag] != NULL)
        {
                free(fragment_list[i_frag]);
		i_frag++;
        }

        delete[] fragment_list;
}

void load_fragments_binary(char* mapped_reads_fp, vector<t_fragment*>* fore_strand_frags, vector<t_fragment*>* rev_strand_frags)
{
	printf("Loading fragments from %s.\n", mapped_reads_fp);

        FILE* f_mapped_reads = fopen(mapped_reads_fp, "rb");
        if(f_mapped_reads == NULL)
        {
                printf("Could not open mapped reads file %s @ %s(%d).\n", mapped_reads_fp, __FILE__, __LINE__);
                exit(0);
        }

	// File reading loop.
	while(1)
	{
		// Read the current fragment index.
		int cur_fragment_base_index;
		if(fread(&cur_fragment_base_index, sizeof(int), 1, f_mapped_reads) != 1)
		{
			break;
		}
		
		// Read the current sequenced fragment length.
		int cur_fragment_sequenced_fragment_length;
                if(fread(&cur_fragment_sequenced_fragment_length, sizeof(int), 1, f_mapped_reads) != 1)
		{
			printf("Binary fragment file %s ended while reading fragment length @ %s(%d)\n", mapped_reads_fp, __FILE__, __LINE__);
			exit(0);
		}

		// Read the direction.
		char cur_strand_char;
		if(fread(&cur_strand_char, sizeof(char), 1, f_mapped_reads) != 1)
		{
                        printf("Binary fragment file %s ended while reading fragment direction @ %s(%d)\n", mapped_reads_fp, __FILE__, __LINE__);
                        exit(0);
		}
	
		//printf("Adding the fragment %d\n", cur_fragment_base_index);
	
		// Allocate the new fragment and add it to the list.
		if(cur_strand_char == 'F')
		{
			// Allocate a new fragment and add it to the forward strand fragment list.
			t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
			new_fragment->base_index = cur_fragment_base_index;
			new_fragment->strand_char = cur_strand_char;
			//new_fragment->sequenced_fragment_length = strlen(effective_fragment);
			new_fragment->sequenced_fragment_length = cur_fragment_sequenced_fragment_length;

			//chr_fragments->push_back(new_fragment);
			fore_strand_frags->push_back(new_fragment);

		} // Forward strand check.
		else if(cur_strand_char == 'R')
		{
			// For fragments on reverse strand, the base index is set to the 5' side of the fragment, where the indexing is with respect to forward strand.
			// Allocate and initialize a fragment and add it to the reverse strand fragment list.
			t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
			//new_fragment->base_index = chr_index + (strlen(cur_fragment) - 1); // Corrected index, note negative 1 correction.
			new_fragment->base_index = cur_fragment_base_index + cur_fragment_sequenced_fragment_length; // Corrected index, note negative 1 correction. [PeakSeq compatilibity]
			new_fragment->strand_char = cur_strand_char;
			//new_fragment->sequenced_fragment_length = strlen(effective_fragment);
			new_fragment->sequenced_fragment_length = cur_fragment_sequenced_fragment_length;
			//chr_fragments->push_back(new_fragment);
			rev_strand_frags->push_back(new_fragment);

		} // Reverse strand check
		else
		{
			
		}
	}

	fclose(f_mapped_reads);

	printf("Loaded %d fragments on forward strand.\n", fore_strand_frags->size());
        printf("Loaded %d fragments on reverse strand.\n", rev_strand_frags->size());

        // First must make sure the overly emnriched sies are removed from the fragment list: Remove the positions for which the count is greater than 3.
        // Basic rule is that the starting indices of tags on forward and reverse strands should not be mixed until enrichment profiles are generated. For
        // example, do not sort the fragments.

        // First get rid of the extra fragments on forward strand.
        vector<t_fragment*>* pruned_fore_fragments = prune_fragment_list_by_strand(fore_strand_frags, 'F', 30000000);
        vector<t_fragment*>* pruned_rev_fragments = prune_fragment_list_by_strand(rev_strand_frags, 'R', 300000000);

        //delete_fragments(fore_strand_frags);
        //delete_fragments(rev_strand_frags);

        for(int i_frag = 0; i_frag < fore_strand_frags->size(); i_frag++)
        {
                free(fore_strand_frags->at(i_frag));
        }
        fore_strand_frags->clear();

        for(int i_frag = 0; i_frag < rev_strand_frags->size(); i_frag++)
        {
                free(rev_strand_frags->at(i_frag));
        }
        rev_strand_frags->clear();

        // Combine all the fragments.
        vector<t_fragment*>* pruned_fragments = new vector<t_fragment*>();
        for(int i_frag = 0; i_frag < pruned_fore_fragments->size(); i_frag++)
        {
/*
                t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
                new_fragment_node->base_index = pruned_fore_fragments->at(i_frag)->base_index;
                new_fragment_node->sequenced_fragment_length = pruned_fore_fragments->at(i_frag)->sequenced_fragment_length;
                new_fragment_node->strand_char = pruned_fore_fragments->at(i_frag)->strand_char;

                fore_strand_frags->push_back(new_fragment_node);
*/
                fore_strand_frags->push_back(pruned_fore_fragments->at(i_frag));
        } // i_fore_frag loop.

        for(int i_frag = 0; i_frag < pruned_rev_fragments->size(); i_frag++)
        {
/*
                t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
                new_fragment_node->base_index = pruned_rev_fragments->at(i_frag)->base_index;
                new_fragment_node->sequenced_fragment_length = pruned_rev_fragments->at(i_frag)->sequenced_fragment_length;
                new_fragment_node->strand_char = pruned_rev_fragments->at(i_frag)->strand_char;

                rev_strand_frags->push_back(new_fragment_node);
*/
                rev_strand_frags->push_back(pruned_rev_fragments->at(i_frag));
        } // i_fore_frag loop.


/*
        for(int i_frag = 0; i_frag < fore_strand_frags->size(); i_frag++)
        {
                printf("[bin fore] %d, %d, %c\n", fore_strand_frags->at(i_frag)->base_index, fore_strand_frags->at(i_frag)->sequenced_fragment_length, fore_strand_frags->at(i_frag)->strand_char);
        }


        for(int i_frag = 0; i_frag < rev_strand_frags->size(); i_frag++)
        {
                printf("[bin rev] %d, %d, %c\n", rev_strand_frags->at(i_frag)->base_index, rev_strand_frags->at(i_frag)->sequenced_fragment_length, rev_strand_frags->at(i_frag)->strand_char);
        }
*/
        // Free fragment memory.
        pruned_fore_fragments->clear();
        delete(pruned_fore_fragments);

        pruned_rev_fragments->clear();
        delete(pruned_rev_fragments);

}

void dump_fragments_binary(char* mapped_reads_bin_fp, vector<t_fragment*>* fragments)
{
	FILE* f_mapped_reads_bin = open_f(mapped_reads_bin_fp, "wb");
	for(int i_frag = 0; i_frag < fragments->size(); i_frag++)
	{
		// It is very important to do the following correction to the fragments on the negative strand: The indices in the file correspond to the left side index
		// of the tag, without regard to it being on the reverse strand. In the loaded memory, the starting index of a fragment points to it 5' end.
		if(fragments->at(i_frag)->strand_char == 'R')
		{
			int rev_frag_base_i = fragments->at(i_frag)->base_index - fragments->at(i_frag)->sequenced_fragment_length;
			fwrite(&rev_frag_base_i, sizeof(int), 1, f_mapped_reads_bin);
		}
		else
		{
			fwrite(&fragments->at(i_frag)->base_index, sizeof(int), 1, f_mapped_reads_bin);
		}
		fwrite(&fragments->at(i_frag)->sequenced_fragment_length, sizeof(int), 1, f_mapped_reads_bin);
		fwrite(&fragments->at(i_frag)->strand_char, sizeof(char), 1, f_mapped_reads_bin);
	} // i_frag loop.
	fclose(f_mapped_reads_bin);
}

void load_fragments(t_chromosome* chr, char* mapped_reads_fp, vector<t_fragment*>* fore_strand_frags, vector<t_fragment*>* rev_strand_frags)
{
	printf("Loading fragments from %s.\n", mapped_reads_fp);

	FILE* f_mapped_reads = fopen(mapped_reads_fp, "r");
	if(f_mapped_reads == NULL)
	{
		printf("Could not open mapped reads file %s @ %s(%d).\n", mapped_reads_fp, __FILE__, __LINE__);
		exit(0);
	}

	char cur_fragment[100];
	char quality_str[20];
	char strand_char;
	int chr_index;

	// Read and validate the mapped reads in the file.
	while(fscanf(f_mapped_reads, "%s %s %c %d", cur_fragment, quality_str, &strand_char, &chr_index) == 4)
	{
		//printf("Read %s\n", cur_fragment);
		char effective_fragment[200];
		memset(effective_fragment, 0, 200 * sizeof(char));
		// Analyze the fragment: Check the leading and following 'N's. This affects the length of the fragment.

		// If a chromosome is specified, do sanity check for the sequence.
		if(chr != NULL)
		{
			// Get rid of the 'N's at the beginning.
			int first_i = 0;
			int last_i = 0;
			for(int i = 0; i < strlen(cur_fragment); i++)
			{
				if(cur_fragment[i] == 'I')
				{
					printf("Found an 'I' in fragment data:\n%s\n", cur_fragment);
					getc(stdin);
				}
	
				if(cur_fragment[i] != 'N' && cur_fragment[i] != 'I')
				{
					first_i = i;
					break;
				}
			} // i loop.

			if(first_i != 0)
			{
				printf("first_i is not 0, i.e., there are N's at the beginning, they will be ignored in a brain dead manner for\n%s\n", cur_fragment);
				//getc(stdin);
			}

			//printf("Current fragment length: %d (%s)\n", strlen(cur_fragment), cur_fragment);	
			for(int i = strlen(cur_fragment)-1; i > 0; i--)
			{
				if(cur_fragment[i] == 'I')
				{
					printf("Found an 'I' in fragment data:\n%s\n", cur_fragment);
					getc(stdin);
				}
	
				if(cur_fragment[i] != 'N' && cur_fragment[i] != 'I')
				{
					last_i = i;
					break;
				}
			} // i loop.

			if(last_i !=  strlen(cur_fragment)-1)
			{
				printf("last i is not strlen(cur_fragment)-1, %d, %d\n", last_i, strlen(cur_fragment)-1);
				//exit(0);
			}
	
			// Copy effective fragment.
			for(int i = first_i; i <= last_i; i++)
			{
				effective_fragment[i - first_i] = cur_fragment[i];
			}
	
			//printf("%s -> %s (%d)         \r", cur_fragment, effective_fragment, strlen(effective_fragment));
	/*
			if(strlen(effective_fragment) != strlen(cur_fragment))
			{
				printf("\n");
			}
	*/
			// Do a sanity check: Compare the subsequence at the indicated index with the fragment.
			char chr_subseq[100];
			memset(chr_subseq, 0, 100 * sizeof(char));

			if(strand_char == 'F')
			{
				for(int i = 0; i < strlen(effective_fragment); i++)
				{
					chr_subseq[i] = chr->x_nuc(chr_index + i);
				} // i loop
		
                                // Note that base index for chr_subseq is 0, although base index for chromosome is 1.
				if(!check_fragment_quality(effective_fragment, chr_subseq, quality_str))
				{
					printf("%s (%c):\n%s\n(%s)\n%s\n", quality_str, strand_char, effective_fragment, cur_fragment,chr_subseq);
					//getc(stdin);
				}

				// Allocate a new fragment and add it to the forward strand fragment list.
	                        t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
	                        new_fragment->base_index = chr_index;
	                        new_fragment->strand_char = strand_char;
	                        //new_fragment->sequenced_fragment_length = strlen(effective_fragment);
	                        new_fragment->sequenced_fragment_length = strlen(cur_fragment);
	
	                        //chr_fragments->push_back(new_fragment);
	                        fore_strand_frags->push_back(new_fragment);

			} // Forward strand check.
			else if(strand_char == 'R')
			{
				memset(chr_subseq, 0, 100 * sizeof(char));
	
				/*
				chr_index-1: Chromosomes nucleotide data is 0-base indexed in t_chromosome object and 1-base indexed in fragment data file.
				strlen(effective_fragment)-1: The highest index that can be used for a string is (strlen(string) - 1).
				*/
				for(int i = 0; i < strlen(effective_fragment); i++)
				{
					char reverse_strand_nuc = chr->x_nuc(chr_index + (strlen(effective_fragment)-1) - i);
					//printf("%c", reverse_strand_nuc);
		
					chr_subseq[i] = get_complement_nuc_by_nuc(reverse_strand_nuc);
				} // i loop
		
				//printf("\n");

				// For fragments on reverse strand, the base index is set to the 5' side of the fragment, where the indexing is with respect to forward strand.
				// Note that base index for chr_subseq is 0, although base index for chromosome is 1.
				if(!check_fragment_quality(effective_fragment, chr_subseq, quality_str))
				{
					printf("%s (%c):\n%s\n%s\n", quality_str, strand_char, effective_fragment, chr_subseq);
					//getc(stdin);
				}
				
				if(first_i  != 0)
				{
					//printf("Matched a case where first_i not 0 and quality ok on reverse strand:\n%s\n%s\n", effective_fragment, chr_subseq);
					//getc(stdin);
				}

				// Allocate and initialize a fragment and add it to the reverse strand fragment list.
	                        t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
	                        //new_fragment->base_index = chr_index + (strlen(cur_fragment) - 1); // Corrected index, note negative 1 correction.
				new_fragment->base_index = chr_index + strlen(cur_fragment); // Corrected index, note negative 1 correction. [PeakSeq compatilibity]
	                        new_fragment->strand_char = strand_char;
	                        //new_fragment->sequenced_fragment_length = strlen(effective_fragment);
	                        new_fragment->sequenced_fragment_length = strlen(cur_fragment);
	
	                        //chr_fragments->push_back(new_fragment);
	                        rev_strand_frags->push_back(new_fragment);
			} // Reverse strand check
		} // Check if the chromosome is specified.
		else
		{
                        if(strand_char == 'F')
                        {
		                t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
		                new_fragment->base_index = chr_index;
		                new_fragment->strand_char = strand_char;
				new_fragment->sequenced_fragment_length = strlen(cur_fragment);
		
		                //chr_fragments->push_back(new_fragment);
				// Adding all of the fragments to forward strand.
	                        fore_strand_frags->push_back(new_fragment);
			}
			else if(strand_char == 'R')
			{
                                //printf("[txt] %d, %d, %c\n", chr_index, strlen(cur_fragment), strand_char);

                                // Allocate and initialize a fragment and add it to the reverse strand fragment list.
                                t_fragment* new_fragment = (t_fragment*)malloc(sizeof(t_fragment));
                                //new_fragment->base_index = chr_index + (strlen(cur_fragment) - 1); // Corrected index, note negative 1 correction.
                                new_fragment->base_index = chr_index + strlen(cur_fragment); // Corrected index, note negative 1 correction. [PeakSeq compatilibity]
                                new_fragment->strand_char = strand_char;
                                //new_fragment->sequenced_fragment_length = strlen(effective_fragment);
                                new_fragment->sequenced_fragment_length = strlen(cur_fragment);

                                rev_strand_frags->push_back(new_fragment);

			}
		} // Chromosome check.
	} // curent fragment data reading loop.

	fclose(f_mapped_reads);

	printf("Loaded %d fragments on forward strand.\n", fore_strand_frags->size());
        printf("Loaded %d fragments on reverse strand.\n", rev_strand_frags->size());

	// First must make sure the overly emnriched sies are removed from the fragment list: Remove the positions for which the count is greater than 3.
	// Basic rule is that the starting indices of tags on forward and reverse strands should not be mixed until enrichment profiles are generated. For 
	// example, do not sort the fragments.
	
	// First get rid of the extra fragments on forward strand.
        vector<t_fragment*>* pruned_fore_fragments = prune_fragment_list_by_strand(fore_strand_frags, 'F', 30000000);
        vector<t_fragment*>* pruned_rev_fragments = prune_fragment_list_by_strand(rev_strand_frags, 'R', 30000000);

	//delete_fragments(fore_strand_frags);
	//delete_fragments(rev_strand_frags);

	for(int i_frag = 0; i_frag < fore_strand_frags->size(); i_frag++)
	{
		free(fore_strand_frags->at(i_frag));
	}
	fore_strand_frags->clear();

        for(int i_frag = 0; i_frag < rev_strand_frags->size(); i_frag++)
        {
                free(rev_strand_frags->at(i_frag));
        }
	rev_strand_frags->clear();

	// Combine all the fragments.
        vector<t_fragment*>* pruned_fragments = new vector<t_fragment*>();
	for(int i_frag = 0; i_frag < pruned_fore_fragments->size(); i_frag++)
	{
/*
		t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
		new_fragment_node->base_index = pruned_fore_fragments->at(i_frag)->base_index;
		new_fragment_node->sequenced_fragment_length = pruned_fore_fragments->at(i_frag)->sequenced_fragment_length;
		new_fragment_node->strand_char = pruned_fore_fragments->at(i_frag)->strand_char;
	
		fore_strand_frags->push_back(new_fragment_node);
*/
                fore_strand_frags->push_back(pruned_fore_fragments->at(i_frag));
	} // i_fore_frag loop.	

        for(int i_frag = 0; i_frag < pruned_rev_fragments->size(); i_frag++)
        {
/*
                t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
                new_fragment_node->base_index = pruned_rev_fragments->at(i_frag)->base_index;
                new_fragment_node->sequenced_fragment_length = pruned_rev_fragments->at(i_frag)->sequenced_fragment_length;
                new_fragment_node->strand_char = pruned_rev_fragments->at(i_frag)->strand_char;

                rev_strand_frags->push_back(new_fragment_node);
*/
		rev_strand_frags->push_back(pruned_rev_fragments->at(i_frag));
        } // i_fore_frag loop.

/*
        for(int i_frag = 0; i_frag < fore_strand_frags->size(); i_frag++)
        {
                printf("[txt fore] %d, %d, %c\n", fore_strand_frags->at(i_frag)->base_index, fore_strand_frags->at(i_frag)->sequenced_fragment_length, fore_strand_frags->at(i_frag)->strand_char);
        }


        for(int i_frag = 0; i_frag < rev_strand_frags->size(); i_frag++)
        {
                printf("[txt rev] %d, %d, %c\n", rev_strand_frags->at(i_frag)->base_index, rev_strand_frags->at(i_frag)->sequenced_fragment_length, rev_strand_frags->at(i_frag)->strand_char);
        }
*/

        // Free fragment memory.
        pruned_fore_fragments->clear();
	delete(pruned_fore_fragments);

	pruned_rev_fragments->clear();
	delete(pruned_rev_fragments);
}

bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str)
{
        int mismatches = 0;
        for(int i = 0; i < strlen(fragment); i++)
        {
                char upper_subseq = toupper(chr_subseq[i]);
                char upper_frag = toupper(fragment[i]);

                if((upper_subseq != upper_frag) &&
                        (upper_subseq != 'N') &&
                        (upper_frag != 'N'))
                {
                        mismatches++;
                        //printf("Fragments are not matching for %s at %d (%s). Enriched fragment is:\n%s\n%s\n", chr_fps->at(i_chr), chr_index, quality_str, cur_seq, chr_subseq);
                        //getc(stdin);
                        //exit(0);
                }
        } // i loop.

        if(mismatches == 0)
        {
		return true;
/*
                if(strcmp(quality_str, "U0") == 0)
                {
                        return(true);
                }
                else
                {
                        return(false);
                }
*/
        }
        else if(mismatches == 1)
        {
                if(strcmp(quality_str, "U1") == 0)
                {
                        return(true);
                }
                else
                {
                        return(false);
                }
        }
        else if(mismatches == 2)
        {
                if(strcmp(quality_str, "U2") == 0)
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

char get_complement_nuc_by_nuc(char nuc_2_complement)
{
        if(nuc_2_complement == 'A' || nuc_2_complement == 'a')
        {
                return('T');
        }
        else if(nuc_2_complement == 'C' || nuc_2_complement == 'c')
        {
                return('G');
        }
        else if(nuc_2_complement == 'G' || nuc_2_complement == 'g')
        {
                return('C');
        }
        else if(nuc_2_complement == 'T' || nuc_2_complement == 't')
        {
                return('A');
        }
        else
        {
                printf("Cannot complement %c\n", nuc_2_complement);
                return(0);
        }
}

bool sort_fragments(t_fragment* frag1, t_fragment* frag2)
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

vector<t_fragment*>* prune_fragment_list_by_strand(vector<t_fragment*>* chr_fragments, char strand_2_prune, int n_max_reps)
{
	// Sort the fragment list, cannpt prune a non-sorted fragment list.
	sort(chr_fragments->begin(), chr_fragments->end(), sort_fragments);

	printf("Pruning %d fragments on '%c' strand with maximum reps of %d.\n", chr_fragments->size(), strand_2_prune, n_max_reps);

        // First get rid of the extra fragments on forward strand.
        int* rep_cnts = new int[chr_fragments->size() + 2];
	memset(rep_cnts, 0, sizeof(int) * (chr_fragments->size() + 1));
        for(int i_frag = 0; i_frag < chr_fragments->size(); i_frag++)
        {
                if(chr_fragments->at(i_frag)->strand_char == strand_2_prune)
                {
                        if(i_frag > 0 &&
				chr_fragments->at(i_frag)->base_index == chr_fragments->at(i_frag-1)->base_index)
                        {
                                rep_cnts[i_frag] = rep_cnts[i_frag-1] + 1;
                        }
                        else // This is a new fragment set its copy number to 1.
                        {
                                rep_cnts[i_frag] = 1;
                        }
                }
		else
		{
			printf("WTF???\n");
			exit(0);
			rep_cnts[i_frag] = 0;
		}
        } // i_frag loop.

        vector<t_fragment*>* pruned_fragments = new vector<t_fragment*>();

        for(int i_frag = 0; i_frag < chr_fragments->size(); i_frag++)
        {
                if(chr_fragments->at(i_frag)->strand_char == strand_2_prune)
                {
	                if(rep_cnts[i_frag] <= n_max_reps)
	                {
	                        t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
	                        new_fragment_node->base_index = chr_fragments->at(i_frag)->base_index;
	                        new_fragment_node->sequenced_fragment_length = chr_fragments->at(i_frag)->sequenced_fragment_length;
	                        new_fragment_node->strand_char = chr_fragments->at(i_frag)->strand_char;
	
	                        pruned_fragments->push_back(new_fragment_node);
        	        }
		}
        } // i_frag loop.

	delete [] rep_cnts;

	printf("Pruned to %d fragments.\n", pruned_fragments->size());

	return(pruned_fragments);
}

/*
Forwardiz'ing: The libraryies store the 5' end of the fragments in the memory, by convention. Forward'izing puts all the bases to 5' side and reverses the fragment
orientation, in addition, it can do the tag extension (only place in the whole library that tag extension is done) and also sorting of the new fragments. This function
returns a new set of fragments that are all on the forward strand with the proper extensions applied, and sorted.
*/
vector<t_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_fragment*>* fore_frag_list, vector<t_fragment*>* rev_frag_list, int enrichment_fragment_length)
{
	vector<t_fragment*>* combined_frags = new vector<t_fragment*>();

	if(fore_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < fore_frag_list->size(); i_frag++)
		{
			t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));
			new_fragment_node->base_index = fore_frag_list->at(i_frag)->base_index;

			// If the enrichment length is greater than the actual seq1uence tag length, update the sequenced length.
			if(enrichment_fragment_length > fore_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				new_fragment_node->sequenced_fragment_length = enrichment_fragment_length;
			}
			else
			{
				new_fragment_node->sequenced_fragment_length = fore_frag_list->at(i_frag)->sequenced_fragment_length;
			}


			new_fragment_node->strand_char = fore_frag_list->at(i_frag)->strand_char;

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	if(rev_frag_list != NULL)
	{
		for(int i_frag = 0; i_frag < rev_frag_list->size(); i_frag++)
		{
			// If this enrichment fragment length is set to 0, this means that there is no fragment extension.

			t_fragment* new_fragment_node = (t_fragment*)malloc(sizeof(t_fragment));

			if(enrichment_fragment_length < rev_frag_list->at(i_frag)->sequenced_fragment_length)
			{
				new_fragment_node->base_index = (rev_frag_list->at(i_frag)->base_index > rev_frag_list->at(i_frag)->sequenced_fragment_length)?(rev_frag_list->at(i_frag)->base_index - rev_frag_list->at(i_frag)->sequenced_fragment_length):(1); 
				new_fragment_node->sequenced_fragment_length = rev_frag_list->at(i_frag)->sequenced_fragment_length;
			}
			else
			{
				new_fragment_node->base_index = (rev_frag_list->at(i_frag)->base_index > enrichment_fragment_length)?(rev_frag_list->at(i_frag)->base_index - enrichment_fragment_length):(1); // Correct the fragment indices for fragment on the reverse strand.
				new_fragment_node->sequenced_fragment_length = enrichment_fragment_length;
			}

			new_fragment_node->strand_char = 'F'; // Change the place of these to forward strand.

			combined_frags->push_back(new_fragment_node);
		} // i_frag loop.
	}

	// Sort one last time.
	sort(combined_frags->begin(), combined_frags->end(), sort_fragments);

	return(combined_frags);
}

int* get_n_reads_per_window(int n_wins, vector<t_fragment*>* frag_list)
{
	int* n_reads_per_window = new int[n_wins + 2];
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		n_reads_per_window[i_win] = 0;
	} // i_win loop

	for(int i_frag = 0; i_frag < frag_list->size(); i_frag++)
	{
		int cur_i_win = frag_list->at(i_frag)->base_index / MEG_BASE;

		n_reads_per_window[cur_i_win]++;
	} // i_nuc loop.

	return(n_reads_per_window);
}

void delete_frag_ends(vector<t_frag_end*>* frag_ends)
{
	for(int i_end = 0; i_end < frag_ends->size(); i_end++)
	{
		delete(frag_ends->at(i_end));
	} // i_end loop.

	delete(frag_ends);
}

void add_fragment_ends(vector<t_frag_end*>* frag_ends, int start, int length)
{
	t_frag_end* frag_start = new t_frag_end();
	frag_start->pos = start;
	frag_start->side = FRAG_START;

	t_frag_end* frag_end = new t_frag_end();
	frag_end->pos = start + length; // Note that this position points to the nucleotide where profile height goes down. Therefore this is 1 nucleotide after the end of the fragment. The fragment ends at (start + length - 1).
	frag_end->side = FRAG_END;

	// Add these ends.
	frag_ends->push_back(frag_start);
	frag_ends->push_back(frag_end);
}

// Return the fragment ends for a valid set of fragments. Note that the tag length extension is not applied here any more.
vector<t_frag_end*>* get_frag_ends(vector<t_fragment*>* fragments)
{
        // Initialize counts.
        vector<t_frag_end*>* frag_ends = new vector<t_frag_end*>();

		//printf("Building profile with %d fragments.\n", fragments->size());

        // Proces each fragment in the fragments list.
        for(int i_frag = 0; i_frag < fragments->size(); i_frag++)
        {
                t_fragment* cur_frag = fragments->at(i_frag);

                if(cur_frag->strand_char == 'F')
                {
					// Add two points: One for beginning of the fragment and one for the end:
					//this->add_fragment_ends(cur_frag->base_index - 1, enrichment_fragment_length);
					add_fragment_ends(frag_ends, cur_frag->base_index, cur_frag->sequenced_fragment_length);
                }
                else if(cur_frag->strand_char == 'R')
                {
					add_fragment_ends(frag_ends, cur_frag->base_index - (cur_frag->sequenced_fragment_length - 1), cur_frag->sequenced_fragment_length);
                }
                else
                {
                        printf("Could not resolve strand char '%c' for current fragment, i_frag = %d(%d) @ %s(%d)\n", cur_frag->strand_char, i_frag, fragments->size(), __FILE__, __LINE__);
                        exit(0);
                }
        } // i_frag loop.

        // Sort the fragment ends.
        sort(frag_ends->begin(), frag_ends->end(), sort_ends);
        //fprintf(stderr, "Done!\n");

		return(frag_ends);
}

bool sort_ends(t_frag_end* end1, t_frag_end* end2)
{
	return(end1->pos < end2->pos);
}


