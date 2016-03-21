#include <stdio.h>
#include <stdlib.h>
#include "annot_region_tools.h"
#include <algorithm>
#include <string.h>
#include <math.h>
#include "rng.h"
#include "seed_manager.h"
#include "utils.h"
#include "ansi_string.h"
#include "nomenclature.h"
#include "genomics_coords.h"

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

bool __DUMP_ANNOT_REGION_TOOLS_MSGS__ = false;

/*
Make sure that the external data structure pointers are correctly set for merging/intersecting operations between the newly allocated regions. This can be ensured by not reallocating
new regions while restructuring the regions per strands/chromosomes. restructure_annot_regions function may be very useful for this.
*/
void dump_BED(char* bed_fp, vector<t_annot_region*>* annot_regions)
{
	FILE* f_bed = open_f(bed_fp, "w");

	for(int i_reg = 0; i_reg < (int)annot_regions->size(); i_reg++)
	{
		if(annot_regions->at(i_reg)->strand == '-' ||
			annot_regions->at(i_reg)->strand == '+')
		{
			if(annot_regions->at(i_reg)->name == NULL)
			{
				annot_regions->at(i_reg)->name = new char[5];
				strcpy(annot_regions->at(i_reg)->name, ".");
			}

			// Translate the start and end to CODEBASE's start and end.
			fprintf(f_bed, "%s\t%d\t%d\t%s\t.\t%c\n", annot_regions->at(i_reg)->chrom, 
				translate_coord(annot_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(annot_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				annot_regions->at(i_reg)->name,
				annot_regions->at(i_reg)->strand);
		}
		else
		{
			fprintf(f_bed, "%s\t%d\t%d\n", annot_regions->at(i_reg)->chrom, 
				translate_coord(annot_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(annot_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base));
		}
	} // i_ref loop.

	fclose(f_bed);
}

void delete_restructured_annot_regions(t_sorted_annot_region_lists* restructured_region_lists)
{
	for(int i_chr = 0; i_chr < (int)restructured_region_lists->chr_ids->size(); i_chr++)
	{
		delete [] restructured_region_lists->chr_ids->at(i_chr);
		delete restructured_region_lists->neg_strand_regions_per_chrom[i_chr];
		delete restructured_region_lists->pos_strand_regions_per_chrom[i_chr];
		delete restructured_region_lists->regions_per_chrom[i_chr];
	} // i_chr loop.

	delete [] restructured_region_lists->neg_strand_regions_per_chrom;
	delete [] restructured_region_lists->pos_strand_regions_per_chrom;
	delete [] restructured_region_lists->regions_per_chrom;
	delete restructured_region_lists->chr_ids;

	delete restructured_region_lists;
}

vector<t_annot_region*>* load_BED(char* bed_fp)
{
	vector<t_annot_region*>* bed_regions = new vector<t_annot_region*>();

	FILE* f_bed = open_f(bed_fp, "r");
	
	while(1)
	{
		char* cur_line = getline(f_bed);
		if(cur_line == NULL)
		{
			break;
		}

		bool skip_line = check_line_skip(cur_line);

		if(!skip_line)
		{
			// Read the mandatory fields.
			int i_cur_char = 0;
			char strand_char = '+';
			char region_name[1000];
			strcpy(region_name, ".");
			int region_score = 0;

			char chrom[100];
			memset(chrom, 0, 100);
			if(!t_string::get_next_token(cur_line, chrom, 100, "\t", i_cur_char))
			{
				fprintf(stderr, "Could not read chromosome in %s\n", cur_line);
				exit(0);
			}
			

			int start = 0;
			char start_str[100];
			memset(start_str, 0, 100);
			if(!t_string::get_next_token(cur_line, start_str, 100, "\t", i_cur_char))
			{
				fprintf(stderr, "Could not read start in %s\n", cur_line);
				exit(0);
			}
			else
			{
				start = atoi(start_str);
			}

			int end = 0;
			char end_str[100];
			memset(end_str, 0, 100);
			if(!t_string::get_next_token(cur_line, end_str, 100, "\t", i_cur_char))
			{
				fprintf(stderr, "Could not read end in %s\n", cur_line);
				exit(0);
			}
			else
			{
				end = atoi(end_str);
			}

			
			// Keep on reading and load until line is finished.
			bool line_finished = false;
			char name_str[1000];
			memset(name_str, 0, 1000);
			if(!t_string::get_next_token(cur_line, name_str, 1000, "\t", i_cur_char))
			{
				// Could not get the next token.
				line_finished = true;
			}
			else
			{
				// Region name is already initialized to '.'.
				strcpy(region_name, name_str);
			}

			if(!line_finished)
			{
				char score_str[1000];
				memset(score_str, 0, 1000);
				if(!t_string::get_next_token(cur_line, score_str, 1000, "\t", i_cur_char))
				{
					line_finished = true;
				}
				else
				{
					//strcpy(region_score, score_str);
					region_score = atoi(score_str);
				}
			}

			if(!line_finished)
			{
				char strand_str[10];
				memset(strand_str, 0, 10);
				if(!t_string::get_next_token(cur_line, strand_str, 10, "\t", i_cur_char))
				{
					line_finished = true;
				}
				else
				{
					strand_char = strand_str[0];
				}
			}

			/*if(sscanf(cur_line, "%s %d %d", chrom, &start, &end) != 3)
			{
				printf("Could not read the mandatory fields from BED file line:\n%s\n", cur_line);
				exit(0);
			}*/

			//char strand = '+'; 
			//if(sscanf(cur_line, "%*s %*s %*s %*s %*s %c", &strand) != 1)
			//{
			//	// Could not read the strand from file, set it to '+' by default.
			//	strand = '+';
			//}

			//int region_score = 0;
			//if(sscanf(cur_line, "%*s %*s %*s %*s %d", &region_score) != 1)
			//{
			//	// Could not read the score from the file.
			//}

			//char region_name[1000];
			//if(sscanf(cur_line, "%*s %*s %*s %s", region_name) != 1)
			//{
			//	// Could not read the strand from file, set it to '+' by default.
			//	strcpy(region_name, ".");
			//}

			t_annot_region* new_region = new t_annot_region();

			// Translate the start and end of the BED coordinates to the codebase coordinates.
			//new_region->start = start - BED_START_BASE + CODEBASE_START_BASE;
			//new_region->end = end - BED_END_BASE + CODEBASE_END_BASE;
			new_region->start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
			new_region->end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
			new_region->chrom = new_region->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(new_region->chrom);
			new_region->name = new char[strlen(region_name) + 2];
			strcpy(new_region->name, region_name);
			new_region->strand = strand_char;
			new_region->data = NULL;
			new_region->score = region_score;
			new_region->sort_info = NULL;
		
			bed_regions->push_back(new_region);
		} // line skip check.

		delete [] cur_line;
	} // file reading loop.

	fclose(f_bed);

	//printf("Loaded %d regions from %s\n", bed_regions->size(), bed_fp);

	if(!validate_region_coords(bed_regions))
	{
		fprintf(stderr, "The coordinates are not valid for %s\n", bed_fp);
		exit(0);
	}

	return(bed_regions);
}

bool check_line_skip(char* cur_line)
{
	// Skip empty lines.
	if(strlen(cur_line) == 0)
	{
		return(true);
	}

	// Check comment.
	if(cur_line[0] == '#')
	{
		return(true);
	}
	
	// Check track info line.
	char* first_word = new char[strlen(cur_line) + 2];
	strcpy(first_word, cur_line);
	sscanf(cur_line, "%s", first_word);

	if(strcmp(first_word, "track") == 0)
	{
		delete [] first_word;
		return(true);
	}

	delete [] first_word;

	return(false);
}

t_annot_region* duplicate_region(t_annot_region* region_2_dup)
{
	t_annot_region* dup_region = new t_annot_region();

	if(region_2_dup->chrom != NULL)
	{
		dup_region->chrom = new char[strlen(region_2_dup->chrom) + 2];
		strcpy(dup_region->chrom, region_2_dup->chrom);
	}
	else
	{
		dup_region->chrom = NULL;
	}

	if(region_2_dup->name != NULL)
	{
		dup_region->name = new char[strlen(region_2_dup->name) + 2];
		strcpy(dup_region->name, region_2_dup->name);
	}
	else
	{
		dup_region->name = NULL;
	}

	dup_region->start = region_2_dup->start;
	dup_region->end = region_2_dup->end;
	dup_region->strand = region_2_dup->strand;
	dup_region->n_exons = dup_region->n_exons;
	dup_region->score = dup_region->score;
	dup_region->thick_start = dup_region->thick_start;
	dup_region->thick_end = dup_region->thick_end;

	dup_region->intervals = NULL;

	return(dup_region);
}

vector<char*>* get_chr_ids(vector<t_annot_region*>* annot_regions)
{
	//printf("/*Determining*/ the chromosome list from region list.\n");
	vector<char*>* chr_ids = new vector<char*>();

	for(int i_reg = 0; i_reg < (int)annot_regions->size(); i_reg++)
	{
		bool is_new_id = true;

		// Check if the the current chromosome id is included in the list of chromosome ids already collected, otherwise add it.
		for(int i_chr = 0; 
			is_new_id && i_chr < (int)chr_ids->size(); 
			i_chr++)
		{
			if(strcmp(chr_ids->at(i_chr), annot_regions->at(i_reg)->chrom) == 0)
			{
				is_new_id = false;
			}
		} // i_chr loop

		if(is_new_id)
		{
			char* new_id = new char[strlen(annot_regions->at(i_reg)->chrom) + 2];
			strcpy(new_id, annot_regions->at(i_reg)->chrom);
			chr_ids->push_back(new_id); // Add the new id.
		}
	} // i_reg loop.

	//printf("Determined %d chromosomes from region list.\n", chr_ids->size());
	return(chr_ids);
}

t_annot_region* get_empty_region()
{
	t_annot_region* new_reg = new t_annot_region();
	new_reg->chrom = NULL;
	new_reg->start = 0;
	new_reg->end = 0;
	new_reg->strand = 0;
	new_reg->intervals = 0;
	new_reg->name = 0;
	new_reg->data = 0;
	new_reg->sort_info = NULL;

	return(new_reg);
}

void delete_annot_regions(t_annot_region* region)
{
	if(region != NULL)
	{
		delete [] region->chrom;

		if(region->name != NULL)
		{
			delete [] region->name;
		}

		delete(region);
	}
}

/*
Clean a list of annotated regions.
*/
void delete_annot_regions(vector<t_annot_region*>* region_list)
{
	for(int i_reg = 0; i_reg < (int)region_list->size(); i_reg++)
	{
		if(region_list->at(i_reg) != NULL)
		{
			delete [] region_list->at(i_reg)->chrom;

			if(region_list->at(i_reg)->name != NULL)
			{
				delete [] region_list->at(i_reg)->name;
			}

			delete(region_list->at(i_reg));
		}
	} // i_reg loop.

	region_list->clear();
	delete(region_list);
}


bool validate_region_coords(vector<t_annot_region*>* regions)
{
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		if(regions->at(i_reg)->start > regions->at(i_reg)->end)
		{
			fprintf(stderr, "%s:%d-%d is not valid.\n", regions->at(i_reg)->chrom, regions->at(i_reg)->start, regions->at(i_reg)->end);
			return(false);
		}
	} // i_reg loop.

	return(true);
}


