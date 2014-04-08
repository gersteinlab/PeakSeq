#include "signal_track_tools.h"
#include "../../../lib/database/annotation/annot_region_tools.h"
#include "../../../lib/database/mapped_read/fragment.h"
#include "../../utils/file/utils.h"
#include "../../utils/ansi_string/ansi_string.h"
#include "../../../Genomics/ChIP_Seq/Peak_Selector/peak_region.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

// t_profile_site structure based processing of signal tracks.
/*
The basic idea behind t_profile_site is simple: Each profile site indicates the updated signal height starting
at the site's index till the next profile site. For example:
chr1	1	0
chr1	100	8
chr1	200	9
chr1	250 0

Shows a signal profile that starts from 100 and goes till 250. As a bedgraph file this would be:
chr1	100	200	8
chr1	200	250	9

All the functionality is provided at vector<t_profile_site*>* level. The loading functions load_bedGraph and load_sgr
loads the files into memory and separates the signal profiles with respect to chromosome and store it in t_signal_profile structure. 
The separated signals can then be processed with the basic functions: For each chromosome, extract the regions of interest by buffer_signal_profile
function, then use aggregation/correlation functions for processing the data.

For the region based functions, the last index in the region is NOT included in processing. For example, while buffering the 
signal profiles, the last index in the region of interest is not copied to the buffered profile.
*/

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

/*
Subtract the first profile from the second and return the remainder profile. This can be negative.
*/
vector<t_profile_site*>* subtract_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2)
{
	vector<t_profile_site*>* delta_prof1 = get_delta_profile_per_profile(profile1);

	vector<t_profile_site*>* delta_prof2 = get_delta_profile_per_profile(profile2);

	// Negate the delta_profile for profile2
	amplify_profile(delta_prof2, -1.0);

	// Concatenate the fragment ends and sort them, then generate the profile with the new list of fragment ends.
	delta_prof1->insert(delta_prof1->end(), delta_prof2->begin(), delta_prof2->end());

	// Sort the fragment ends.
	sort(delta_prof1->begin(), delta_prof1->end(), sort_profile_sites);
		
	// Get the profile site for the aggregated fragment ends.
	vector<t_profile_site*>* aggregate_profile = get_profile_per_delta_profile(delta_prof1);

	// Free memory.
	delete_profile_sites(delta_prof1);
	delete(delta_prof2);

	return(aggregate_profile);
}

/*
Compute the variance for the list of profiles.
*/
vector<t_profile_site*>* get_variance_profile(vector<vector<t_profile_site*>*>* profile_list)
{
	if(profile_list->size() == 1)
	{
		return(NULL);
	}

	// Aggregate the profiles then compute the mean.
	vector<t_profile_site*>* mean_profile = aggregate_profiles(profile_list);

	// Get the mean: Divide by total number of profiles.
	amplify_profile(mean_profile, 1.0 / (double)(profile_list->size()));

	// Estimate the sample variance: Subtract the mean profile from each profile, then take the corrected mean for unbiased variance estimation.
	vector<vector<t_profile_site*>*>* squared_sample_mean_shifted_profiles = new vector<vector<t_profile_site*>*>();

	// Subtract the sample mean from each profile and store it.
	for(int i_p = 0; i_p < profile_list->size(); i_p++)
	{
		vector<t_profile_site*>* cur_mean_shifted_profile = subtract_profiles(profile_list->at(i_p), mean_profile);

		// Square the current mean shifted profile.
		square_profile(cur_mean_shifted_profile);
		squared_sample_mean_shifted_profiles->push_back(cur_mean_shifted_profile);
	} // i_p loop.
	
	// Aggregate the squared/mean shifted profiles.
	vector<t_profile_site*>* variance_profile = aggregate_profiles(squared_sample_mean_shifted_profiles);

	// Divide the profile by corrected factor to estimate unbiased variance.
	amplify_profile(variance_profile, 1.0 / (double)(variance_profile->size() - 1));

	// Clean memory.
	for(int i_prof = 0; i_prof < squared_sample_mean_shifted_profiles->size(); i_prof++)
	{
		delete_profile_sites(squared_sample_mean_shifted_profiles->at(i_prof));
	} // i_prof loop.
	delete(squared_sample_mean_shifted_profiles);

	// Free mean profile memory.
	delete_profile_sites(mean_profile);

	return(variance_profile);
}

void square_profile(vector<t_profile_site*>* profile)
{
	for(int i_p = 0; i_p < profile->size(); i_p++)
	{
		double sq_val = profile->at(i_p)->height * profile->at(i_p)->height;
		profile->at(i_p)->height = sq_val;
	}
}

void amplify_profile(vector<t_profile_site*>* profile, double amp)
{
	for(int i_site = 0; i_site < profile->size(); i_site++)
	{
		profile->at(i_site)->height *= amp;
	} // i_site loop.
}

// Take two signal profiles and sum them up: Get the fragment ends and take union, then convert the union to signal profile.
vector<t_profile_site*>* aggregate_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2)
{
	//printf("%d profile1.\n", profile1->size());
	//for(int i = 0; i < profile1->size(); i++)
	//{
	//	printf("%d\t%lf\n", profile1->at(i)->i_nuc, profile1->at(i)->height);
	//}

	//printf("%d profile2.\n", profile2->size());
	//for(int i = 0; i < profile2->size(); i++)
	//{
	//	printf("%d\t%lf\n", profile2->at(i)->i_nuc, profile2->at(i)->height);
	//}

	//printf("Aggreagating profiles (%d, %d).\n", profile1->size(), profile2->size());
	//fprintf(stderr, "Getting delta profile for profile1\n");
	vector<t_profile_site*>* delta_prof1 = get_delta_profile_per_profile(profile1);
	//printf("%d delta_profile1.\n", delta_prof1->size());
	//for(int i = 0; i < delta_prof1->size(); i++)
	//{
	//	printf("%d\t%lf\n", delta_prof1->at(i)->i_nuc, delta_prof1->at(i)->height);
	//}

	//fprintf(stderr, "Getting delta profile for profile2\n");
	vector<t_profile_site*>* delta_prof2 = get_delta_profile_per_profile(profile2);
	//printf("%d delta_profile2.\n", delta_prof2->size());
	//for(int i = 0; i < delta_prof2->size(); i++)
	//{
	//	printf("%d\t%lf\n", delta_prof2->at(i)->i_nuc, delta_prof2->at(i)->height);
	//}

	// Concatenate the fragment ends and sort them, then generate the profile with the new list of fragment ends.
	//fprintf(stderr, "Inserting delta profiles.\n");
	delta_prof1->insert(delta_prof1->end(), delta_prof2->begin(), delta_prof2->end());

	// Sort the fragment ends.
	//fprintf(stderr, "Sorting delta profiles.\n");
	sort(delta_prof1->begin(), delta_prof1->end(), sort_profile_sites);
		
	// Get the profile site for the aggregated fragment ends.
	//fprintf(stderr, "Final aggregation.\n");
	vector<t_profile_site*>* aggregate_profile = get_profile_per_delta_profile(delta_prof1);

	// Free memory.
	delete_profile_sites(delta_prof1);

	// Note that the profile sites for second delta profile is deleted with the delete_profile_sites call above.
	delete(delta_prof2);

	return(aggregate_profile);
}

vector<t_profile_site*>* aggregate_profiles(vector<vector<t_profile_site*>*>* profile_lists)
{
	//printf("Aggreagating profiles (%d, %d).\n", profile1->size(), profile2->size());
	vector<t_profile_site*>* total_delta_prof = new vector<t_profile_site*>();

	for(int i_prof = 0; i_prof < profile_lists->size(); i_prof++)
	{
		vector<t_profile_site*>* cur_delta_prof = get_delta_profile_per_profile(profile_lists->at(i_prof));

		total_delta_prof->insert(total_delta_prof->end(), cur_delta_prof->begin(), cur_delta_prof->end());
	} // i_prof loop.	

	// Sort the fragment ends.
	sort(total_delta_prof->begin(), total_delta_prof->end(), sort_profile_sites);
		
	// Get the profile site for the aggregated fragment ends.
	vector<t_profile_site*>* aggregate_profile = get_profile_per_delta_profile(total_delta_prof);

	// Free memory.
	delete_profile_sites(total_delta_prof);

	return(aggregate_profile);
}


void delete_profile_sites(vector<t_profile_site*>* profile_sites)
{
	for(int i_site = 0; i_site < profile_sites->size(); i_site++)
	{
		delete(profile_sites->at(i_site));
	} // i_site loop.

	delete(profile_sites);
}

/*
Load a wig file and dump the sgr file from that. The wig file is parsed with respect to the chromosome identity.
*/
void parse_WIG_formatted_signal_track(vector<char*>* chr_ids, char* parsed_signal_tracks_op_dir, char* wig_fp)
{
	FILE* f_wig = fopen(wig_fp, "r");
	if(f_wig == NULL)
	{
		printf("Could not open WIG file %s.\n", wig_fp);
		exit(0);
	}

	// These are the buffered profile sites per chromosome. 
	vector<vector<t_wig_entry*>*>* wig_entries_per_chrom = new vector<vector<t_wig_entry*>*>();
    for(int i = 0; i < chr_ids->size(); i++)
    {
		wig_entries_per_chrom->push_back(new vector<t_wig_entry*>());
    }

	// Line buffer for loading the file.
	char cur_line[1000];

	int start = 0;
	int step = 1;
	int span = 1;
	int signal_nuc_start = 0;
	double signal_value;
	char chrom[100];

	while(parse_new_WIG_entry(f_wig,
						wig_entries_per_chrom,
						chr_ids,
						chrom, 
						start,
						step,
						span,
						signal_nuc_start, 
						signal_value))
	{
		// Process next entry in the wig file.
	}

	// Sort the wig file entries for each chromosome. This is necessary in case that the tracks in the wig file are put in random order. The order in the sgr file should be sorted correctly.
    printf("Sorting the wig entries per chromosome.\n");
	for(int i = 0; i < chr_ids->size(); i++)
    {
		sort(wig_entries_per_chrom->at(i)->begin(), wig_entries_per_chrom->at(i)->end(), sort_wig_entries);

		printf("%d entries in %s\n", wig_entries_per_chrom->at(i)->size(), chr_ids->at(i));
	}

	// Open the files.
    vector<FILE*>* signal_track_f_ptrs = new vector<FILE*>();
    for(int i = 0; i < chr_ids->size(); i++)
    {		
		// Check if there are entries for this chromosome.
		if(wig_entries_per_chrom->at(i)->size() > 0)
		{
	        char new_fn[100];

			sprintf(new_fn, "%s/%s_signal_track.sgr", parsed_signal_tracks_op_dir, chr_ids->at(i));
				signal_track_f_ptrs->push_back(open_f(new_fn, "w"));
		}
		else
		{
			signal_track_f_ptrs->push_back(NULL);
		}
    } // i loop for chromosomes


	// Write the files for each chromosome, whose entry list is formed above.
	for(int i = 0; i < chr_ids->size(); i++)
    {
		if(wig_entries_per_chrom->at(i)->size() > 0)
		{
			// For this chromosome, convert the wig entries into sgr profile sites then dump the sgr entries.
			dump_sgr_per_WIG_entries(wig_entries_per_chrom->at(i), signal_track_f_ptrs->at(i), chr_ids->at(i));
		}
	} // i loop for chromosomes

    for(int i = 0; i < chr_ids->size(); i++)
    {
		if(signal_track_f_ptrs->at(i) != NULL)
		{
			fclose(signal_track_f_ptrs->at(i));
		}
    } // i loop.

	fclose(f_wig);
}

void dump_sgr_per_WIG_entries(vector<t_wig_entry*>* wig_entries, FILE* f_sgr, char* chr_id)
{
	/* 
	Go over all the wig entries and generate the sgr profile.
	The idea is to start from the first index in the wig entries and update the sgr entries per signal profile change.
	*/
	vector<t_profile_site*>* sgr_prof_sites = new vector<t_profile_site*>();

	double cur_height = 0; // This is the current height.
	int cur_wig_pos = 0;
	for(int i_wig_site = 0; i_wig_site < wig_entries->size(); i_wig_site++)
	{
		t_wig_entry* cur_wig_entry = wig_entries->at(i_wig_site);

		// Process all the height changes between the cur_pos and the new wig_site.
		if(cur_wig_entry->pos == cur_wig_pos + 1)
		{
			// The wig entry is following cur_pos, if there is no hiehgt change, ignore this.
			if(cur_wig_entry->value != cur_height)
			{
				t_profile_site* new_sgr_prof_site = new t_profile_site();
				new_sgr_prof_site->height = cur_wig_entry->value;
				new_sgr_prof_site->i_nuc = cur_wig_entry->pos; // Do not forget the index base relocations for SGR and WIG files.

				sgr_prof_sites->push_back(new_sgr_prof_site);
				fprintf(f_sgr, "%s\t%d\t%lf\n", chr_id, new_sgr_prof_site->i_nuc, new_sgr_prof_site->height);
			}
		}
		else if(cur_wig_entry->pos > cur_wig_pos + 1)
		{
			//printf("Adding a 0 height: Jumped from %d to %d\n", cur_wig_pos, cur_wig_entry->pos);

			// The wig entry is further away from the last processed position: Add a 0 height first. This is necessary for ensuring that sgr file shows the decrease of height to 0. 
			t_profile_site* new_sgr_prof_site = new t_profile_site();
			new_sgr_prof_site->height = 0.0; // The height drops to 0.
			new_sgr_prof_site->i_nuc = cur_wig_pos + 1; // Do not forget the index base relocations for SGR and WIG files.

			sgr_prof_sites->push_back(new_sgr_prof_site);
			fprintf(f_sgr, "%s\t%d\t%lf\n", chr_id, new_sgr_prof_site->i_nuc, new_sgr_prof_site->height);

			// Then Add the next entry as the new height.
			new_sgr_prof_site = new t_profile_site();
			new_sgr_prof_site->height = cur_wig_entry->value;
			new_sgr_prof_site->i_nuc = cur_wig_entry->pos; // Do not forget the index base relocations for SGR and WIG files.

			sgr_prof_sites->push_back(new_sgr_prof_site);
			fprintf(f_sgr, "%s\t%d\t%lf\n", chr_id, new_sgr_prof_site->i_nuc, new_sgr_prof_site->height);
		}
		else
		{
			printf("WTF is going on? cur_pos: %d, new_wig_entry_pos: %d\n", cur_wig_pos, cur_wig_entry->pos);
			exit(0);
		}

		// Update the current height and current wig position.
		cur_wig_pos = cur_wig_entry->pos + cur_wig_entry->span - 1;
		cur_height = cur_wig_entry->value;
	} // i_wig_site loop.
}

vector<t_wig_entry*>* get_sig_track_list_pointer_by_chr_file_name(vector<vector<t_wig_entry>*>* sig_track_lists, vector<char*>* chr_ids, char* chr_fn)
{
	return(NULL);
}

vector<t_wig_entry*>* get_sig_track_list_pointer_by_chr_id(vector<vector<t_wig_entry*>*>* sig_track_lists, vector<char*>* chr_ids, char* chr_id)
{
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		if(t_string::compare_strings_ci(chr_ids->at(i_chr), chr_id))
		{
			return(sig_track_lists->at(i_chr));
		}
	} // i_chr loop.

	printf("Could not find the wig entry list for chromosome with id %s", chr_id);
	exit(0);
	return(NULL);
}

FILE* get_sig_track_file_pointer_by_chr_file_name(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_fn)
{
	return(NULL);
}

FILE* get_sig_track_file_pointer_by_chr_id(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_id)
{
	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		if(t_string::compare_strings_ci(chr_ids->at(i_chr), chr_id))
		{
			return(fragment_f_ptrs->at(i_chr));
		}
	} // i_chr loop.

	printf("Could not find the wig entry list for chromosome with id %s", chr_ids);
	exit(0);
	return(NULL);
}

bool sort_wig_entries(t_wig_entry* entry1, t_wig_entry* entry2)
{
	return(entry1->pos < entry2->pos);
}

bool sort_bedGraph_entries(t_bedGraph_entry* entry1, t_bedGraph_entry* entry2)
{
	return(entry1->start < entry2->start);
}

// This function updates all the chromosome, start, span, step information. Return false when there is no more data to read from the file.
bool parse_new_WIG_entry(FILE* f_wig,
						vector<vector<t_wig_entry*>*>* sig_track_lists,
						vector<char*>* chr_ids,
						char* chrom, 
						int& start,
						int& step,
						int& span,
						int& signal_nuc_start, 
						double& signal_value)
{
	char cur_line[1000];

	// Load the file: Assuming that there is one track in the file.
	// The first line is track name: Load the track name.
	if(x_fgets(cur_line, 1000, f_wig) == NULL)
	{
		printf("File ended.\n");
		return(false);
	}

	// Read and skip the header lines.
	int l_line = strlen(cur_line);
	while(cur_line[l_line-1] == '\\' || cur_line[0] == '#')
	{
		// Read the next line.
		if(x_fgets(cur_line, 1000, f_wig) == NULL)
		{
			printf("The WIG file ended while reading a track line or after a comment.\n");
			exit(0);
		}

		l_line = strlen(cur_line);
	} // Header line loading loop.

	char cur_line_id[100];
	sscanf(cur_line, "%s", cur_line_id);

	// This is either a data line, or a new track info line, check if the first string is step declaration.
	if(strcmp(cur_line_id, "track") == 0)
	{
		// A new track is started.
		// Do not process this line unless it is necessary to keep the track information.
		printf("Found track information: %s\n", cur_line);
		return(true);
	}
	else if(strcmp(cur_line_id, "fixedStep") == 0)
	{
		// Read the remaining parameters: step size, span, ...
		// Parse the line with respect to space, read the parameters.
		printf("WIG file contains a track with fixedStep. Not implemented, yet.\n");
		exit(0);
	}
	else if(strcmp(cur_line_id, "variableStep") == 0)
	{
		// Update the chromosome, step, span information, note that these are utilized while buffering the profile site lists.
		parse_WIG_info_line(cur_line, 
							chrom,
							start, 
							step, 
							span);

		printf("Found variable stepping section: %s\n", cur_line);

		// Note that this function will be called again.
		return(true);
	}
	else
	{
		// This is a signal line for a locus, parse the signal value and the locus.
		int _signal_nuc_start;
		double _signal_value;
		if(sscanf(cur_line, "%d %lf", &_signal_nuc_start, &_signal_value) != 2)
		{
			printf("Could not parse the signal line: %s\n", cur_line);
		}
		
		// Copy the values.
		signal_nuc_start = _signal_nuc_start;
		signal_value = _signal_value;

		// Get the correct chromosome and add the wig entry.
		vector<t_wig_entry*>* cur_chr_wig_site_list = get_sig_track_list_pointer_by_chr_id(sig_track_lists, chr_ids, chrom);

		t_wig_entry* new_site = new t_wig_entry();
		new_site->pos = signal_nuc_start;
		new_site->value = _signal_value;
		new_site->span = span; // This is kinda important while processing this entry.
		cur_chr_wig_site_list->push_back(new_site);

		//printf("New WIG entry: %s@%d-%d:%lf\n", chrom, signal_nuc_start, signal_nuc_start + span - 1, new_site->value);

		return(true);
	}
}

void parse_WIG_info_line(char* wig_info_line, 
						char* chrom,
						int& start, 
						int& step, 
						int& span)
{
	start = 1;
	step = 1;
	span = 1;

	t_string* line_str = new t_string(wig_info_line);

	t_string_tokens* info_line_tokens = line_str->tokenize_by_chars(" =");

	for(int i_tok = 0; i_tok < info_line_tokens->size(); i_tok++)
	{
		if(info_line_tokens->at(i_tok)->compare_ci("chrom"))
		{
			i_tok++; // Update token index.
			strcpy(chrom, info_line_tokens->at(i_tok)->str());
		}
		else if(info_line_tokens->at(i_tok)->compare_ci("start"))
		{
			i_tok++;
			start = atoi(info_line_tokens->at(i_tok)->str());
		}
		else if(info_line_tokens->at(i_tok)->compare_ci("step"))
		{
			i_tok++;
			step = atoi(info_line_tokens->at(i_tok)->str());
		}
		else if(info_line_tokens->at(i_tok)->compare_ci("span"))
		{
			i_tok++;
			span = atoi(info_line_tokens->at(i_tok)->str());
		}
	} // i_tok loop.

	line_str->clean_tokens(info_line_tokens);
	delete(line_str);

	printf("WIG file: %s, starting at %d with steps of %d, and span of %d.\n", chrom, start, step, span);
}


void parse_bedGraph_formatted_signal_track(vector<char*>* chr_ids, char* parsed_signal_tracks_op_dir, char* bedgraph_fp)
{
	FILE* f_bedgraph = fopen(bedgraph_fp, "r");
	if(f_bedgraph == NULL)
	{
		printf("Could not open bedGraph file %s.\n", bedgraph_fp);
		exit(0);
	}

	// These are the buffered profile sites per chromosome. 
	vector<vector<t_wig_entry*>*>* wig_entries_per_chrom = new vector<vector<t_wig_entry*>*>();
    for(int i = 0; i < chr_ids->size(); i++)
    {
		wig_entries_per_chrom->push_back(new vector<t_wig_entry*>());
    }

	// Line buffer for loading the file.
	char cur_line[1000];

	int start = 0;
	int step = 1;
	int span = 1;
	int signal_nuc_start = 0;
	double signal_value;
	char chrom[100];

	int i_line = 0;
	while(parse_new_bedGraph_entry(f_bedgraph,
									wig_entries_per_chrom,
									chr_ids,
									chrom, 
									start,
									step,
									span,
									signal_nuc_start, 
									signal_value))
	{
		// Process next entry in the wig file.
		i_line++;
		if(i_line % (1000 * 1000) == 0)
		{
			fprintf(stderr, "Processing %d. entry                       \r", i_line);
		}
	}

	// Sort the wig file entries for each chromosome. This is necessary in case that the tracks in the wig file are put in random order. The order in the sgr file should be sorted correctly.
    printf("Sorting the wig entries per chromosome.\n");
	for(int i = 0; i < chr_ids->size(); i++)
    {
		sort(wig_entries_per_chrom->at(i)->begin(), wig_entries_per_chrom->at(i)->end(), sort_wig_entries);

		printf("%d entries in %s\n", wig_entries_per_chrom->at(i)->size(), chr_ids->at(i));
	}

	// Open the files.
    vector<FILE*>* signal_track_f_ptrs = new vector<FILE*>();
    for(int i = 0; i < chr_ids->size(); i++)
    {		
		// Check if there are entries for this chromosome.
		if(wig_entries_per_chrom->at(i)->size() > 0)
		{
	        char new_fn[100];

			sprintf(new_fn, "%s/%s_signal_track.sgr", parsed_signal_tracks_op_dir, chr_ids->at(i));
				signal_track_f_ptrs->push_back(open_f(new_fn, "w"));
		}
		else
		{
			signal_track_f_ptrs->push_back(NULL);
		}
    } // i loop for chromosomes


	// Write the files for each chromosome, whose entry list is formed above.
	for(int i = 0; i < chr_ids->size(); i++)
    {
		if(wig_entries_per_chrom->at(i)->size() > 0)
		{
			// For this chromosome, convert the wig entries into sgr profile sites then dump the sgr entries.
			dump_sgr_per_WIG_entries(wig_entries_per_chrom->at(i), signal_track_f_ptrs->at(i), chr_ids->at(i));
		}
	} // i loop for chromosomes

    for(int i = 0; i < chr_ids->size(); i++)
    {
		if(signal_track_f_ptrs->at(i) != NULL)
		{
			fclose(signal_track_f_ptrs->at(i));
		}
    } // i loop.

	fclose(f_bedgraph);
}


// This function updates all the chromosome, start, span, step information. Return false when there is no more data to read from the file.
bool parse_new_bedGraph_entry(FILE* f_bedgraph,
								vector<vector<t_wig_entry*>*>* sig_track_lists,
								vector<char*>* chr_ids,
								char* chrom, 
								int& start,
								int& step,
								int& span,
								int& signal_nuc_start, 
								double& signal_value)
{
	char cur_line[1000];

	// Load the file: Assuming that there is one track in the file.
	// The first line is track name: Load the track name.
	if(x_fgets(cur_line, 1000, f_bedgraph) == NULL)
	{
		printf("File ended.\n");
		return(false);
	}

	// Read and skip the header lines.
	int l_line = strlen(cur_line);
	while(cur_line[l_line-1] == '\\' || cur_line[0] == '#')
	{
		// Read the next line.
		if(x_fgets(cur_line, 1000, f_bedgraph) == NULL)
		{
			printf("The WIG file ended while reading a track line or after a comment.\n");
			exit(0);
		}

		l_line = strlen(cur_line);
	} // Header line loading loop.

	char cur_line_id[100];
	sscanf(cur_line, "%s", cur_line_id);

	// This is either a data line, or a new track info line, check if the first string is step declaration.
	if(strcmp(cur_line_id, "track") == 0)
	{
		// A new track is started.
		// Do not process this line unless it is necessary to keep the track information.
		printf("Found track information: %s\n", cur_line);
		return(true);
	}
	else
	{
		// This is a signal line for a locus, parse the signal value and the locus.
		int start;
		int end;
		double signal_value;
		char chrom[100];
		if(sscanf(cur_line, "%s %d %d %lf", chrom, &start, &end, &signal_value) != 4)
		{
			printf("Could not parse the signal line: %s\n", cur_line);
			exit(0);
		}
		
		// Get the correct chromosome and add the wig entry.
		if(signal_value > 0.0)
		{
			vector<t_wig_entry*>* cur_chr_wig_site_list = get_sig_track_list_pointer_by_chr_id(sig_track_lists, chr_ids, chrom);

			t_wig_entry* new_site = new t_wig_entry();
			new_site->pos = start;
			new_site->value = signal_value;
			new_site->span = end - start; // end is NOT included in the range; the signal_value is the signal value for [start, end-1].
			cur_chr_wig_site_list->push_back(new_site);

			//printf("New WIG entry: %s@%d-%d:%lf\n", chrom, signal_nuc_start, signal_nuc_start + span - 1, new_site->value);
		}

		return(true);
	}
}

vector<char*>* get_bedGraph_chr_ids(char* bedGraph_fp)
{
	FILE* f_bgr = open_f(bedGraph_fp, "r");
	
	vector<char*>* chr_ids = new vector<char*>();

	char cur_line[1000];
	while(x_fgets(cur_line, 1000, f_bgr))
	{
		char cur_chrom[100];
		if(sscanf(cur_line, "%s %*d %*d %*s", cur_chrom) != 1)
		{
			printf("Could not read the chromosome id from the line %s\n", cur_line);
			exit(0);
		}

		// Check if this id is already added, otw add it.
		bool found_id = false;
		for(int i_id = 0; 
			!found_id && i_id < chr_ids->size(); 
			i_id++)
		{
			if(t_string::compare_strings_ci(chr_ids->at(i_id), cur_chrom))
			{
				found_id = true;
			}
		} // i_id loop.

		if(!found_id)
		{
			char* new_id = new char[strlen(cur_chrom) + 2];
			strcpy(new_id, cur_chrom);
			chr_ids->push_back(new_id);
		}
	} // file reading loop.

	fclose(f_bgr);

	return(chr_ids);
}

t_signal_profile* allocate_signal_profile()
{
	t_signal_profile* new_signal_profile = new t_signal_profile();
	new_signal_profile->chr_ids = new vector<char*>();
	new_signal_profile->profiles_per_chr = new vector<vector<t_profile_site*>*>();

	return(new_signal_profile);
}

vector<char*>* get_sgr_chr_ids(char* sgr_fp)
{
	FILE* f_sgr = open_f(sgr_fp, "r");
	
	vector<char*>* chr_ids = new vector<char*>();

	char cur_line[1000];
	while(x_fgets(cur_line, 1000, f_sgr))
	{
		char cur_chrom[100];
		if(sscanf(cur_line, "%s %*d %*d", cur_chrom) != 1)
		{
			printf("Could not read the chromosome id from the line %s\n", cur_line);
			exit(0);
		}

		// Check if this id is already added, otw add it.
		bool found_id = false;
		for(int i_id = 0; 
			!found_id && i_id < chr_ids->size(); 
			i_id++)
		{
			if(t_string::compare_strings_ci(chr_ids->at(i_id), cur_chrom))
			{
				found_id = true;
			}
		} // i_id loop.

		if(!found_id)
		{
			char* new_id = new char[strlen(cur_chrom) + 2];
			strcpy(new_id, cur_chrom);
			chr_ids->push_back(new_id);
		}
	} // file reading loop.

	fclose(f_sgr);

	return(chr_ids);
}

//void dump_avg_signal_per_windows(vector<char*>* chr_ids, 
//								t_signal_profile_list* profile_sites_per_chrom, 
//								char* fp, 
//								int win_size)
//{
//	//printf("Loaded %d bedGraph entries.\n", bg_entries->size());
//	FILE* f_avg_signal = fopen(fp, "w");
//
//	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
//	{
//		printf("Dumping average signal for %s\n", chr_ids->at(i_chr));
//
//		vector<t_profile_site*>* profile_sites = profile_sites_per_chrom->at(i_chr);
//
//		// Process all the entries accumulated so far.
//		int i_cur_win = 0; // The index for current window being processed.
//		double cur_win_sum = 0.0;
//		int last_block_end = 0; // This is the end of the last block. Note that the next block must start from last_block_end+1 or higher.
//		double cur_height = 0.0;
//
//		//double win_size_in_kb = (double)win_size / 1000;
//
//		// Also the loops must set 0 values for windows that were skipped.
//		for(int i_entry = 0; i_entry < profile_sites->size(); i_entry++)
//		{
//			// Process the signal value for the block: The update is for the 
//			cur_height = (i_entry > 0)?(profile_sites->at(i_entry-1)->height):0.0;
//			for(int i_nuc = last_block_end+1; i_nuc < profile_sites->at(i_entry)->i_nuc; i_nuc++)
//			{
//				int i_win = i_nuc / win_size;
//
//				//printf("%d(%d)\n", i_nuc, i_win);
//
//				// Are we still in the last window? Then just check for maximum update.
//				if(i_win == i_cur_win)
//				{
//					// Update the total height in this window.
//					cur_win_sum += cur_height;
//				}
//				else if(i_win == i_cur_win+1)
//				{
//					// Dump the current maximum value and reset the maximum value for next window: Take the average by the length of the window.
//					double cur_win_avg = cur_win_sum / win_size;
//
//					// Dump only the non-zero entries, saves lots of disk space.
//					if(cur_win_avg > 0.000001)
//					{
//						//fprintf(f_avg_signal, "%s\t%d\t%d\t%.6f\n", chr_ids->at(i_chr), i_cur_win * win_size, ((i_cur_win+1) * win_size), cur_win_avg);
//						fprintf(f_avg_signal, "%.3f\n", cur_win_avg);
//					}
//					else
//					{
//						fprintf(f_avg_signal, "0\n");
//					}
//
//					// Update window counter.
//					i_cur_win++;
//
//					// Initialize the height for next window.
//					cur_win_sum = cur_height;
//				}
//				else
//				{
//					printf("The window index increased unexpectedly: %d to %d for %d to %d\n", i_cur_win, i_win, last_block_end+1, profile_sites->at(i_entry)->i_nuc-1);
//					exit(0);
//				}
//			} // process the nucleotides between the last nucleotide and the start of this block, there are set to 0.
//
//			// Set the last block end for processing next block.
//			last_block_end = profile_sites->at(i_entry)->i_nuc-1;
//		} // i_entry loop.
//	} // i_chr_loop
//
//	fclose(f_avg_signal);
//}

/*
Following is useful for correlating only the entries of the vectors such that at least one component has a value.
*/
double correlate_non_zeros(vector<double>* bgr1_sigs_per_bin, vector<double>* bgr2_sigs_per_bin)
{
	if(bgr1_sigs_per_bin->size() != bgr2_sigs_per_bin->size())
	{
		printf("Cannot correlate the signals with different lengths (%d, %d) @ %s(%d)\n", bgr1_sigs_per_bin->size(), bgr2_sigs_per_bin->size(), __FILE__, __LINE__);
		exit(0);
	}

	double cur_corr = 0.0;

	vector<double>* integrated_signal1_values_per_genome = new vector<double>();
	vector<double>* integrated_signal2_values_per_genome = new vector<double>();

	// Pick the values that are larger than 0.
	for(int i_bin = 0; i_bin < bgr1_sigs_per_bin->size(); i_bin++)
	{
		if(bgr1_sigs_per_bin->at(i_bin) > 0.001 || bgr2_sigs_per_bin->at(i_bin) > 0.001)
		{
			integrated_signal1_values_per_genome->push_back(bgr1_sigs_per_bin->at(i_bin));
			integrated_signal2_values_per_genome->push_back(bgr2_sigs_per_bin->at(i_bin));
		}
	} // i_chr loop.

	double corr = correlate(integrated_signal1_values_per_genome, integrated_signal2_values_per_genome);

	delete(integrated_signal1_values_per_genome);
	delete(integrated_signal2_values_per_genome);

	return(corr);
}

void binarize_bedGraph(char* fp, char* op_fp)
{
	t_signal_profile* signal_profile = load_bedGraph(fp);

	dump_binary_profile(signal_profile, op_fp);

	// Free memory.
	delete_signal_profile(signal_profile);
}

void dump_binary_profile(t_signal_profile* signal_profile, 
						char* op_fp)
{
	FILE* f_op_fp = open_f(op_fp, "wb");

	// For each chromosome in the list, dump the chromosome id, followed by all the profile sites.
	for(int i_chr = 0; i_chr < signal_profile->chr_ids->size(); i_chr++)
	{
		vector<t_profile_site*>* cur_chr_prof = signal_profile->profiles_per_chr->at(i_chr);
		// Dump the chromosome id.

		// Dump a 10 byte for each chromosome id.
		char chr_entry[10];
		memset(chr_entry, 0, 10);
		strcpy(chr_entry, signal_profile->chr_ids->at(i_chr));
		fwrite(chr_entry, sizeof(char), 10, f_op_fp);

		for(int i_site = 0; i_site < cur_chr_prof->size(); i_site++)
		{
			fwrite(&(cur_chr_prof->at(i_site)->i_nuc), sizeof(int), 1, f_op_fp);
			fwrite(&(cur_chr_prof->at(i_site)->height), sizeof(double), 1, f_op_fp);
		} // i_site loop.

		// End the entries for this chromosome.
		int i_end = 0;
		double end_height = -2.0;

		fwrite(&(i_end), sizeof(int), 1, f_op_fp);
		fwrite(&(end_height), sizeof(double), 1, f_op_fp);
	} // i_chr loop.

	fclose(f_op_fp);
}

t_signal_profile* load_binary_profile(char* fp)
{
	t_signal_profile* signal_profile = new t_signal_profile();

	signal_profile->chr_ids = new vector<char*>();
	signal_profile->profiles_per_chr = new t_signal_profile_list();

	FILE* f_fp = open_f(fp, "rb");
	while(1)
	{
		char _new_chr[11];
		if(fread(_new_chr, sizeof(char), 10, f_fp) != 10)
		{
			break;
		}
		else
		{
			printf("Adding %s\n", _new_chr);
		}

		// Add the new chromosome.
		char* new_chr = new char[strlen(_new_chr) + 1];
		strcpy(new_chr, _new_chr);
		signal_profile->chr_ids->push_back(new_chr);

		// Add a new list.
		vector<t_profile_site*>* cur_chr_profile_sites = new vector<t_profile_site*>();

		// Read the entries for this chromosome.
		while(1)
		{
			int cur_pos;
			if(fread(&cur_pos, sizeof(int), 1, f_fp) != 1)
			{
				printf("Could not read the position.");
				exit(0);
			}

			double cur_height;
			if(fread(&cur_height, sizeof(double), 1, f_fp) != 1)
			{
				printf("Could not read the height.");
				exit(0);
			}

			// Found the end of the entries for this chromosome.
			if(cur_height < 0)
			{
				break;
			}
			else
			{
				t_profile_site* new_site = new t_profile_site();
				new_site->i_nuc = cur_pos;
				new_site->height = cur_height;

				// Add the new site to the last list.
				cur_chr_profile_sites->push_back(new_site);
			}

			//printf("%d\t%lf\n", cur_pos, cur_height);
		} // site loading loop for this chromosome.

		// Add the new list to the signal profile.
		signal_profile->profiles_per_chr->push_back(cur_chr_profile_sites);
		printf("%d sites in %s\n", cur_chr_profile_sites->size(), new_chr);
	} // file reading loop.

	fclose(f_fp);

	return(signal_profile);
}

/*
chr_ids is a restricted set the set of chromosome ids that is different from the chromosome ids of the signal profile.
*/
void dump_avg_signal_per_windows(vector<char*>* chr_ids,
								t_signal_profile* signal_profile,
								char* fp, 
								int win_size)
{
	int corr_start = 0;
	int corr_end = 250 * 1000 * 1000;

	// Set the intervals for computing the correlations for each profile: Generate the regions that will be used to bin the profiles.
	vector<t_annot_region*>* bins = new vector<t_annot_region*>();
	int cur_start = corr_start;
	while((cur_start + win_size) < corr_end)
	{
		t_annot_region* new_region = new t_annot_region();
		new_region->start = cur_start;
		new_region->end = cur_start + win_size; // Note that this is going to be excluded in the extracted signal profile.
		new_region->chrom = new char[10];
		strcpy(new_region->chrom, "chr1");
		bins->push_back(new_region);
		cur_start += win_size;
	} // bin generation.

	//printf("Loaded %d bedGraph entries.\n", bg_entries->size());
	FILE* f_avg_signal = fopen(fp, "w");

	for(int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		printf("Dumping average signal for %s\n", chr_ids->at(i_chr));

		int i_sig_chr = get_i_chr(signal_profile->chr_ids, chr_ids->at(i_chr));

		if(i_sig_chr < signal_profile->chr_ids->size())
		{
			//vector<double>* integrated_signal_per_bins = get_integrated_signals_per_bins(cur_profile, i_min, i_max, win_size);
			vector<double>* integrated_signals_per_bins = new vector<double>();

			// Buffer the profiles for the list of windows that are generated for both profiles.
			t_signal_profile_list* bins_profiles = buffer_profile_per_list(signal_profile->profiles_per_chr->at(i_sig_chr), bins, 1.0);

			// Compute the integration of the signals over the windows
			for(int i_win = 0; i_win < bins_profiles->size(); i_win++)
			{
				integrated_signals_per_bins->push_back(integrate_signal(bins_profiles->at(i_win)));
				//printf("%d-%d: %lf\n", profile1->at(0)->i_nuc + (i_win * bin_size), profile1->at(0)->i_nuc + ((i_win + 1)* bin_size), profile1_signals_per_wins->back());
			} // i_win loop.

			for(int i_win = 0; i_win < integrated_signals_per_bins->size(); i_win++)
			{
				if(integrated_signals_per_bins->at(i_win) > 0.001)
				{
					fprintf(f_avg_signal, "%.3f\n", integrated_signals_per_bins->at(i_win));
				}
				else
				{
					fprintf(f_avg_signal, "0\n", integrated_signals_per_bins->at(i_win));
				}
			} // i_win loop.

			// Free memory.
			for(int i_win = 0; i_win < bins_profiles->size(); i_win++)
			{
				delete_profile_sites(bins_profiles->at(i_win));
			}
			delete(bins_profiles);

			delete(integrated_signals_per_bins);
		}
		else
		{
			for(int i_win = 0; i_win < bins->size(); i_win++)
			{
				fprintf(f_avg_signal, "0\n");
			} // i_win loop.
		}
	} // i_chr_loop

	// Delete bins.
	delete_annot_regions(bins);

	fclose(f_avg_signal);
}

void dump_profile(t_signal_profile* signal_profile,
					char* fp,
					int min_i_nuc, 
					int max_i_nuc, 
					double scaling_factor)
{
	FILE* f_prof = open_f(fp, "w");
	if(f_prof == NULL)
	{
		printf("Could not open %s for writing profile data.\n", fp);
		exit(0);
	}

	for(int i_chr = 0; i_chr < signal_profile->profiles_per_chr->size(); i_chr++)
	{
		double height = 0;
		int prev_pos = 0;
		double prev_height = 0;
		int cur_pos = 1;
		vector<t_profile_site*>* profile_sites = signal_profile->profiles_per_chr->at(i_chr);
		for(int i_end = 0; i_end < profile_sites->size(); i_end++)
		{
			/*
			A new region is added between prev_pos and this->frag_ends->at(i_end). If the requested region overlaps
			with the added region, dump the previous position. Note that the height in the added region is the hiehgt
			that is updated at the previous position.
			*/
			// Update the current height: Note that there can be multiple ends at this position.
			height = profile_sites->at(i_end)->height; // The height is scaled below.
			cur_pos = profile_sites->at(i_end)->i_nuc;

			if(cur_pos >= max_i_nuc || 
				i_end == profile_sites->size()-1)
			{
			// The profile passed the requested region completely, can return from the function, dump the previous position and height and return, this is necessary to account for the last position.
				fprintf(f_prof, "%s\t%d\t%lf\n", signal_profile->chr_ids->at(i_chr), cur_pos, (scaling_factor * (double)height));

				break;
			}
			else if(cur_pos < min_i_nuc)
			{
				// Still did not reach the requested region.
			}
			else if(cur_pos >= min_i_nuc)
			{
				fprintf(f_prof, "%s\t%d\t%lf\n", signal_profile->chr_ids->at(i_chr), cur_pos, (scaling_factor * (double)height));
			}	
		} // i_profile site loop.
	} // i_chr loop.
	fclose(f_prof);
}

void dump_profile(vector<t_profile_site*>* profile_sites,
										char* fp, 
										char* chr_id,
										int min_i_nuc, 
										int max_i_nuc, 
										double scaling_factor)
{
	double height = 0;
	FILE* f_prof = open_f(fp, "w");
	if(f_prof == NULL)
	{
		printf("Could not open %s for writing profile data.\n", fp);
		exit(0);
	}

	int cur_pos = 1;
    for(int i_end = 0; i_end < profile_sites->size(); i_end++)
    {
		/*
		A new region is added between prev_pos and this->frag_ends->at(i_end). If the requested region overlaps
		with the added region, dump the previous position. Note that the height in the added region is the hiehgt
		that is updated at the previous position.
		*/
        // Update the current height: Note that there can be multiple ends at this position.
		height = profile_sites->at(i_end)->height; // The height is scaled below.
		cur_pos = profile_sites->at(i_end)->i_nuc;

		if(cur_pos >= max_i_nuc || 
			i_end == profile_sites->size()-1)
		{
			// The profile passed the requested region completely, can return from the function, dump the previous position and height and return, this is necessary to account for the last position.
			fprintf(f_prof, "%s\t%d\t%lf\n", chr_id, cur_pos, (scaling_factor * (double)height));
			break;
		}
		else if(cur_pos < min_i_nuc)
		{
			// Still did not reach the requested region.
		}
		else if(cur_pos >= min_i_nuc)
		{
			// The newly added region overlaps with the requested region, dump the previous region and height, since 
			// the region corresponds to the height of the revious pos.
			fprintf(f_prof, "%s\t%d\t%lf\n", chr_id, cur_pos, (scaling_factor * (double)height));
		}
    } // i_frag
	fclose(f_prof);
}

vector<double>* get_integrated_signals_per_bins(vector<t_profile_site*>* sites, t_annot_region* region, int bin_size)
{
	vector<double>* integrated_signals_per_bins = new vector<double>();

	// Set the intervals for computing the correlations for each profile: Generate the regions that will be used to bin the profiles.
	vector<t_annot_region*>* bins = new vector<t_annot_region*>();
	int cur_start = region->start;
	while((cur_start + bin_size) <= region->end)
	{
		t_annot_region* new_region = new t_annot_region();
		new_region->start = cur_start;
		new_region->end = cur_start + bin_size; // Note that this is going to be excluded in the extracted signal profile.
		new_region->strand = region->strand; // Strand is important since negative strand signals are buffered in the opposite direction compared to positive strand signals.
		new_region->chrom = new char[strlen(region->chrom) + 2];
		strcpy(new_region->chrom, region->chrom);
		bins->push_back(new_region);
		cur_start += bin_size;
	} // bin generation.

	// Buffer the profiles for the list of windows that are generated for both profiles.
	t_signal_profile_list* bins_profiles = buffer_profile_per_list(sites, bins, 1.0);

	// Compute the integration of the signals over the windows
	for(int i_win = 0; i_win < bins_profiles->size(); i_win++)
	{
		integrated_signals_per_bins->push_back(integrate_signal(bins_profiles->at(i_win)));
		//printf("%d-%d: %lf\n", profile1->at(0)->i_nuc + (i_win * bin_size), profile1->at(0)->i_nuc + ((i_win + 1)* bin_size), profile1_signals_per_wins->back());
	} // i_win loop.

	// It is necessary to revert the integrated signals in case the binned region is on negative strand.
	if(region->strand == '-')
	{
		reverse(integrated_signals_per_bins->begin(), integrated_signals_per_bins->end());
	}

	// Free memory.
	for(int i_win = 0; i_win < bins_profiles->size(); i_win++)
	{
		delete_profile_sites(bins_profiles->at(i_win));
	}
	delete(bins_profiles);

	delete_annot_regions(bins);

	return(integrated_signals_per_bins);
}

/*
Correlate two profiles over the win_size bins.
*/
double correlate_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2, int bin_size)
{
	if(profile1->size() == 0 && profile2->size() == 0)
	{
		return(0.0);
	}
	else if(profile1->size() == 0 || profile2->size() == 0)
	{
		printf("Cannot correlate profiles of different width.\n");
		exit(0);
	}

	if((profile1->back()->i_nuc - profile1->at(0)->i_nuc) != (profile2->back()->i_nuc - profile2->at(0)->i_nuc))
	{
		printf("Cannot correlate profiles of different width.\n");
		exit(0);
	}

	//vector<double>* profile1_signals_per_wins = new vector<double>();
	//vector<double>* profile2_signals_per_wins = new vector<double>();

	double correlation = 0.0;

	int corr_start = MAX(profile1->at(0)->i_nuc, profile2->at(0)->i_nuc);
	int corr_end = MIN(profile1->back()->i_nuc, profile2->back()->i_nuc);
	t_annot_region* corr_reg = new t_annot_region();
	corr_reg->chrom = new char[10];
	strcpy(corr_reg->chrom, "chr1");
	corr_reg->strand = '+';
	corr_reg->start = corr_start;
	corr_reg->end = corr_end;

	vector<double>* profile1_signals_per_wins = get_integrated_signals_per_bins(profile1, corr_reg, bin_size);
	vector<double>* profile2_signals_per_wins = get_integrated_signals_per_bins(profile2, corr_reg, bin_size);

	//// Set the intervals for computing the correlations for each profile: Generate the regions that will be used to bin the profiles.
	//vector<t_annot_region*>* bins = new vector<t_annot_region*>();
	//int cur_start = profile1->at(0)->i_nuc;
	//while((cur_start + bin_size) < profile1->back()->i_nuc)
	//{
	//	t_annot_region* new_region = new t_annot_region();
	//	new_region->start = cur_start;
	//	new_region->end = cur_start + bin_size; // Note that this is going to be excluded in the extracted signal profile.
	//	new_region->chrom = new char[10];
	//	strcpy(new_region->chrom, "chr1");
	//	bins->push_back(new_region);
	//	cur_start += bin_size;
	//} // bin generation.

	//// Buffer the profiles for the list of windows that are generated for both profiles.
	//t_signal_profile_list* profile1_wins_profiles = buffer_profile_per_list(profile1, bins, 1.0);
	//t_signal_profile_list* profile2_wins_profiles = buffer_profile_per_list(profile2, bins, 1.0);

	//// Compute the integration of the signals over the windows
	//vector<double>* profile1_signals_per_wins = new vector<double>();
	////printf("Profile 1 signals per window:\n");
	//for(int i_win = 0; i_win < profile1_wins_profiles->size(); i_win++)
	//{
	//	profile1_signals_per_wins->push_back(integrate_signal(profile1_wins_profiles->at(i_win)));
	//	//printf("%d-%d: %lf\n", profile1->at(0)->i_nuc + (i_win * bin_size), profile1->at(0)->i_nuc + ((i_win + 1)* bin_size), profile1_signals_per_wins->back());
	//} // i_win loop.

	//vector<double>* profile2_signals_per_wins = new vector<double>();
	////printf("Profile 2 signals per window:\n");
	//for(int i_win = 0; i_win < profile2_wins_profiles->size(); i_win++)
	//{
	//	profile2_signals_per_wins->push_back(integrate_signal(profile2_wins_profiles->at(i_win)));
	//	//printf("%d-%d: %lf\n", profile2->at(0)->i_nuc + (i_win * bin_size), profile2->at(0)->i_nuc + ((i_win + 1)* bin_size), profile2_signals_per_wins->back());
	//} // i_win loop.

	// Sanity check.
	if(profile1_signals_per_wins->size() != profile2_signals_per_wins->size())
	{
		printf("%s (%d): %d-%d\n", __FILE__, __LINE__, profile1_signals_per_wins->size(), profile2_signals_per_wins->size());
		exit(0);
	}

	// Use the signals per window to compute the correlation: Compute the covariance terms then divide by the estimates of variance.
	correlation = correlate(profile1_signals_per_wins, profile2_signals_per_wins);

	// Free memory: The buffered regions and the window profiles.
	delete(profile1_signals_per_wins);
	delete(profile2_signals_per_wins);

	//for(int i = 0; i < profile1_wins_profiles->size(); i++)
	//{
	//	delete_profile_sites(profile1_wins_profiles->at(i));
	//} // i loop.
	//delete(profile1_wins_profiles);

	//for(int i = 0; i < profile2_wins_profiles->size(); i++)
	//{
	//	delete_profile_sites(profile2_wins_profiles->at(i));
	//} // i loop.
	//delete(profile2_wins_profiles);

	//delete_annot_regions(bins);

	// Return the correlation.
	return(correlation);
}

double correlate_profiles(t_signal_profile* profile1, t_signal_profile* profile2, int bin_size)
{
	double correlation = 0.0;

	vector<double>* profile1_signals_per_wins = new vector<double>();
	vector<double>* profile2_signals_per_wins = new vector<double>();

	for(int i_chr1 = 0; i_chr1 < profile1->chr_ids->size(); i_chr1++)
	{
		int i_chr2 = get_i_chr(profile2->chr_ids, profile1->chr_ids->at(i_chr1));

		if(i_chr2 < profile2->chr_ids->size())
		{
			printf("Correlating signals on %s.\n", profile1->chr_ids->at(i_chr1));

			vector<t_profile_site*>* cur_profile1 = profile1->profiles_per_chr->at(i_chr1);
			vector<t_profile_site*>* cur_profile2 = profile2->profiles_per_chr->at(i_chr2);

			/* 
			Take the smaller of the regions in the chromosomes and buffer it in both profiles: This is necessary since the correlation requires the regions to 
			be of the same length.
			*/
			int corr_start = MAX(cur_profile1->at(0)->i_nuc, cur_profile2->at(0)->i_nuc);
			int corr_end = MIN(cur_profile1->back()->i_nuc, cur_profile2->back()->i_nuc);

			//printf("%s: %d-%d\n", profile1->chr_ids->at(i_chr1), corr_start, corr_end);
			t_annot_region* corr_reg = new t_annot_region();
			corr_reg->chrom = new char[strlen(profile1->chr_ids->at(i_chr1)) + 2];
			strcpy(corr_reg->chrom, profile1->chr_ids->at(i_chr1));
			corr_reg->strand = '+';
			corr_reg->start = corr_start;
			corr_reg->end = corr_end;

			// Compute the integrated signal for each profile over the windows.
			vector<double>* profile1_integrated_signal = get_integrated_signals_per_bins(cur_profile1, corr_reg, bin_size);
			vector<double>* profile2_integrated_signal = get_integrated_signals_per_bins(cur_profile2, corr_reg, bin_size);

			// Add the vectors.
			profile1_signals_per_wins->insert(profile1_signals_per_wins->end(), profile1_integrated_signal->begin(), profile1_integrated_signal->end());
			profile2_signals_per_wins->insert(profile2_signals_per_wins->end(), profile2_integrated_signal->begin(), profile2_integrated_signal->end());

			// Free the vector memories for the integrated signals over bins.
			delete(profile1_integrated_signal);
			delete(profile2_integrated_signal);
		} // chr2 index check.
	} // i_chr1 loop.

	// Sanity check.
	if(profile1_signals_per_wins->size() != profile2_signals_per_wins->size())
	{
		printf("%s (%d): %d-%d\n", __FILE__, __LINE__, profile1_signals_per_wins->size(), profile2_signals_per_wins->size());
		exit(0);
	}

	//printf("Correlating %d windows.\n", profile1_signals_per_wins->size());

	//for(int i = 0; i < profile1_signals_per_wins->size(); i++)
	//{
	//	printf("%d -> %lf\n", i, profile1_signals_per_wins->at(i));
	//} // i loop.

	//for(int i = 0; i < profile1_signals_per_wins->size(); i++)
	//{
	//	printf("%d -> %lf\n", i, profile2_signals_per_wins->at(i));
	//} // i loop.

	// Use the signals per window to compute the correlation: Compute the covariance terms then divide by the estimates of variance.
	correlation = correlate(profile1_signals_per_wins, profile2_signals_per_wins);

	// Free memory: The buffered regions and the window profiles.
	delete(profile1_signals_per_wins);
	delete(profile2_signals_per_wins);

	// Return the correlation.
	return(correlation);
}

double correlate(vector<double>* sig1, vector<double>* sig2)
{
	if(sig1->size() != sig2->size())
	{
		printf("Cannot correlate signals of different length (%d, %d)\n", sig1->size(), sig2->size());
		exit(0);
	}

	double mean1 = 0.0;
	double mean2 = 0.0;

	for(int i = 0; i < sig1->size(); i++)
	{
		mean1 += sig1->at(i);
	} // i loop.
	mean1 /= sig1->size();

	for(int i = 0; i < sig2->size(); i++)
	{
		mean2 += sig2->at(i);
	} // i loop.
	mean2 /= sig2->size();

	double var1 = 0.0;
	double var2 = 0.0;

	for(int i = 0; i < sig1->size(); i++)
	{
		var1 += (sig1->at(i) - mean1) * (sig1->at(i) - mean1);
	} // i loop.
	var1 /= (sig1->size() - 1);
	double std1 = pow(var1, .5);

	for(int i = 0; i < sig2->size(); i++)
	{
		var2 += (sig2->at(i) - mean2) * (sig2->at(i) - mean2);
	} // i loop.
	var2 /= (sig2->size() - 1);
	double std2 = pow(var2, .5);

	double cross_corr = 0.0;
	for(int i = 0; i < sig1->size(); i++)
	{
		cross_corr += (sig1->at(i) - mean1) * (sig2->at(i) - mean2);
	} // i loop.

	// Correction for the sample size in the variance computation.
	cross_corr *= (1 / (double)(sig1->size() - 1));

	// Do a check on the standard deviations.
	if(std1 > 0.000001 && std2 > 0.000001)
	{
		return(cross_corr / (std1 * std2));
	}
	else
	{
		return(0.0);
	}
}

double integrate_signal(vector<t_profile_site*>* profile)
{
	double integral = 0.0;

	for(int i_site = 0; i_site < profile->size()-1; i_site++)
	{
		integral += (profile->at(i_site+1)->i_nuc - profile->at(i_site)->i_nuc) * profile->at(i_site)->height;
	} // i_site loop.

	return(integral);
}

void translate_profiles_to_origin(t_signal_profile_list* profiles)
{
	for(int i_prof = 0; i_prof < profiles->size(); i_prof++)
	{
		vector<t_profile_site*>* cur_profile = profiles->at(i_prof);

		int i_trans = cur_profile->at(0)->i_nuc;
		for(int i_site = 0; i_site < cur_profile->size(); i_site++)
		{
			cur_profile->at(i_site)->i_nuc -= i_trans; // Move the current profile site.
		} // i_site loop.
	} // i_prof loop.
}

/*
After buffering, it is definitely necessary to have height information at start and at end.
*/
vector<t_profile_site*>* buffer_profile_per_region(vector<t_profile_site*>* profile_sites,
													int start, int end,
													double scaling_factor)
{
	//for(int i_ent = 0; i_ent < profile_sites->size(); i_ent++)
	int i_site = 0;

	vector<t_profile_site*>* cur_region_profile = new vector<t_profile_site*>();

	// Move the site pointer forward till next region is found.
	bool localized_site = false;
	while(i_site < profile_sites->size() &&
		start >= profile_sites->at(i_site)->i_nuc)
	{
		// If the next site is after the start of the region, localization is successful.
		if(i_site < profile_sites->size() - 1 &&
			start < profile_sites->at(i_site+1)->i_nuc)
		{
			localized_site++;
			break;
		}
		i_site++;
	} // i_site loop.

	if(localized_site)
	{
		// Normal processing: Add the starting height, move the pointer till the end.
		t_profile_site* new_site = new t_profile_site();
		new_site->height = profile_sites->at(i_site)->height;
		//new_site->i_nuc = profile_sites->at(i_site)->i_nuc;
		new_site->i_nuc = start;
		cur_region_profile->push_back(new_site);
		i_site++;

		while(i_site < profile_sites->size() && 
			end > profile_sites->at(i_site)->i_nuc)
		{
			//fprintf(f_cur_peak_prof, "%s\t%d\t%d\n", annot_regions->at(i_reg)->chrom, profile_sites->at(i_site)->i_nuc, profile_sites->at(i_site)->height);
			t_profile_site* new_site = new t_profile_site();
			new_site->height = profile_sites->at(i_site)->height;
			new_site->i_nuc = profile_sites->at(i_site)->i_nuc;
			cur_region_profile->push_back(new_site);

			i_site++;
		} // loop for processing current peak.

		// Add the 0 value at the end.
		new_site = new t_profile_site();
		new_site->height = 0;
		new_site->i_nuc = end;
		cur_region_profile->push_back(new_site);
	}
	else
	{
		t_profile_site* new_site = new t_profile_site();
		new_site->height = 0;
		new_site->i_nuc = start;
		cur_region_profile->push_back(new_site);
		i_site++;

		// Add the 0 value at the end.
		new_site = new t_profile_site();
		new_site->height = 0;
		new_site->i_nuc = end;
		cur_region_profile->push_back(new_site);
	}

	return(cur_region_profile);
}

/*
Buffers the signal profile (as t_profile_site vectors) for a list of regions.
*/
t_signal_profile_list* buffer_profile_per_list(vector<t_profile_site*>* profile_sites,
															vector<t_annot_region*>* annot_regions,
															double scaling_factor)
{
	t_signal_profile_list* profile_sites_per_regions = new vector<vector<t_profile_site*>*>();

	// Sort the peak annotated regions.
	sort(annot_regions->begin(), annot_regions->end(), sort_regions);
	
	int i_site = 0;

	for(int i_reg = 0; i_reg < annot_regions->size(); i_reg++)
	{
		// At this point, move the i_site pointer back until it moves completely to the left of current region: It may be possible that there is no
		// site that is completely to the left of region, which is still ok since site will not be localized in the following.
		while(i_site < profile_sites->size() &&
			i_site > 0 &&
			annot_regions->at(i_reg)->start <= profile_sites->at(i_site)->i_nuc)
		{
			i_site--;
		} // i_site initialization loop.

		//printf("i = %d\n", i_site);

		// The signal profile list for the current region.
		vector<t_profile_site*>* cur_region_profile = new vector<t_profile_site*>();

		// The aim is to locate i_site that is located inside the region and the height of the profile at the beginning of the region.
		// 
		double start_height = 0.0; // This is the height at the beginning of the region.
		while(i_site < profile_sites->size() &&
			annot_regions->at(i_reg)->start >= profile_sites->at(i_site)->i_nuc)
		{
			// Update the height.
			start_height = profile_sites->at(i_site)->height;

			// Move to the next site, the start_height is still the height just before this region.
			i_site++;
		} // i_site loop.

		//printf("start_height: %lf, i_site: %d.\n", start_height, i_site);

		// At this point, start_height is the height at the beginning of the region, and i_site points to the first profile_site that is to the right of beginning of the region, 
		// which may be outside the region.
		t_profile_site* new_site = new t_profile_site();
		new_site->height = start_height;
		new_site->i_nuc = annot_regions->at(i_reg)->start;
		cur_region_profile->push_back(new_site);

		while(i_site < profile_sites->size() && 
			annot_regions->at(i_reg)->end > profile_sites->at(i_site)->i_nuc)
		{
			//fprintf(f_cur_peak_prof, "%s\t%d\t%d\n", annot_regions->at(i_reg)->chrom, profile_sites->at(i_site)->i_nuc, profile_sites->at(i_site)->height);
			//printf("::%s\t%d\t%lf\n", annot_regions->at(i_reg)->chrom, profile_sites->at(i_site)->i_nuc, profile_sites->at(i_site)->height);
			t_profile_site* new_site = new t_profile_site();
			new_site->height = profile_sites->at(i_site)->height;
			new_site->i_nuc = profile_sites->at(i_site)->i_nuc;
			cur_region_profile->push_back(new_site);

			i_site++;
		} // loop for processing current peak.

		// Add the 0 value at the end.
		new_site = new t_profile_site();
		new_site->height = 0;
		new_site->i_nuc = annot_regions->at(i_reg)->end;
		cur_region_profile->push_back(new_site);

		// Check if the region is on the negative strand, in which case the signal profile must be flipped.
		profile_sites_per_regions->push_back(cur_region_profile);

		// Reverse the profile if the region is on the negative strand.
		if(annot_regions->at(i_reg)->strand == '-')
		{
			//printf("Reverting <%d-%d>\n", annot_regions->at(i_reg)->start, annot_regions->at(i_reg)->end);
			reverse_profile(cur_region_profile);
		}
	} // i_reg loop.

	return(profile_sites_per_regions);
}

void translate_profile(vector<t_profile_site*>* profile, int new_start)
{
	if(profile->size() > 0)
	{
		int delta_i = new_start - profile->at(0)->i_nuc;
		for(int i = 0; i < profile->size(); i++)
		{
			profile->at(i)->i_nuc += delta_i;
		} // i loop.
	}
}

void check_profile(vector<t_profile_site*>* profile_sites)
{
	if(profile_sites->size() == 1)
	{
		printf("The profile is problematic: There is only one site in the profile.\n");
		exit(0);
	}
}

void reverse_profile(vector<t_profile_site*>* profile_sites)
{
	if(profile_sites->size() == 0)
	{
		return;
	}

	if(profile_sites->size() == 1)
	{
		printf("The profile is problematic: There is only one site in the profile.\n");
		exit(0);
	}

	vector<t_profile_site*>* reverted_profile = new vector<t_profile_site*>();

	int i_start = profile_sites->at(0)->i_nuc;
	int i_cur = i_start;

	for(int i = profile_sites->size()-1; i >= 0; i--)
	{
		// i is the index of the element that is going to be processed in this iteration.
		double new_height = 0.0;
		
		// This sets the height for current index.
		if(i > 0)
		{
			new_height = profile_sites->at(i-1)->height;
		}

		t_profile_site* new_site = new t_profile_site();
		new_site->height = new_height;
		new_site->i_nuc = i_cur;
		reverted_profile->push_back(new_site);

		// Following updates the index for the next iteration.
		if(i > 0)
		{
			int new_i = i_cur + (profile_sites->at(i)->i_nuc - profile_sites->at(i-1)->i_nuc);
			//printf("new_i = %d\n", new_i);
			i_cur = new_i;
		}
	} // i loop.

	// Free memory.
	for(int i = 0; i < profile_sites->size(); i++)
	{
		delete(profile_sites->at(i));
	}
	profile_sites->clear();

	// Copy the the reverted profile.
	for(int i = 0; i < reverted_profile->size(); i++)
	{
		profile_sites->push_back(reverted_profile->at(i));
	} // i loop.

	//printf("Reverted profile has %d entries.\n", profile_sites->size());
}

/*
Normalize signal width.
*/
void normalize_signal_profile_length(vector<t_profile_site*>* profile, int total_normalized_width)
{
	sort(profile->begin(), profile->end(), sort_profile_sites);

	int first_i_nuc = profile->at(0)->i_nuc;
	int last_i_nuc = profile->back()->i_nuc;
	double normalization_factor = (double)total_normalized_width / (last_i_nuc - first_i_nuc);

	//printf("%d-%d: %lf\n", first_i_nuc, last_i_nuc, normalization_factor);

	// Translate the profile.
	t_signal_profile_list* profile_list = new t_signal_profile_list();
	profile_list->push_back(profile);
	translate_profiles_to_origin(profile_list);

	for(int i_site = 0; i_site < profile->size(); i_site++)
	{
		double cur_pos = profile->at(i_site)->i_nuc * normalization_factor;
		profile->at(i_site)->i_nuc = (int)cur_pos;
	} // i_site loop.

	// Go over the profile sites and make sure that the sites that overlap with each other are merged (via averaging the heights), this may be the case when
	// a large region is shrinked into a small region.
	vector<t_profile_site*>* new_profile = new vector<t_profile_site*>();
	int i_site = 0; 
	while(i_site < profile->size())
	{
		int cur_pos = profile->at(i_site)->i_nuc;
		double cur_height = 0.0;
		int n_pos = 0;
		while(i_site < profile->size() &&
			profile->at(i_site)->i_nuc == cur_pos)
		{
			//cur_height += profile->at(i_site)->height;
			cur_height = profile->at(i_site)->height;
			i_site++;
			n_pos++;

			//if(n_pos > 1)
			//	printf("Averaging @ %d\n", cur_pos);
		}

		// Take the average of all the heights that are merged.
		//cur_height /= n_pos;

		t_profile_site* new_site = new t_profile_site();
		new_site->height = cur_height;
		new_site->i_nuc = (cur_pos + first_i_nuc); // Translate back to original origin.

		new_profile->push_back(new_site);
	} // i_site loop.

	// Copy the new profile to profile list in the argument.
	for(int i_site = 0; i_site < profile->size(); i_site++)
	{
		delete(profile->at(i_site));
	} // i_site loop.
	profile->clear();

	// Copy the profile sites in the new profile to profile.
	for(int i_site = 0; i_site < new_profile->size(); i_site++)
	{
		profile->push_back(new_profile->at(i_site));
	} 

	// Clean memory.
	new_profile->clear();
	delete(new_profile);
}

void dump_profile_per_list(vector<t_profile_site*>* profile_sites,
												char* dump_dir,
												vector<t_annot_region*>* annot_regions,
												double scaling_factor)
{
	// Sort the peak annotated regions.
	sort(annot_regions->begin(), annot_regions->end(), sort_regions);
	
	//for(int i_ent = 0; i_ent < profile_sites->size(); i_ent++)
	int i_site = 0;

	for(int i_reg = 0; i_reg < annot_regions->size(); i_reg++)
	{
		//// Move the site pointer forward till next region is found.
		//while(i_site < profile_sites->size() &&
		//	annot_regions->at(i_reg)->start > profile_sites->at(i_site)->i_nuc)
		//{
		//	i_site++;
		//}

		// Move the site pointer forward till next region is found.
		while(i_site > 0 &&
			i_site < profile_sites->size() &&
			annot_regions->at(i_reg)->start <= profile_sites->at(i_site)->i_nuc)
		{
			i_site--;
		}

		// Localize the site around the beginning of the region: Find the site that is before start such that the following site is after start.
		bool localized_site = false;
		while(i_site < profile_sites->size() &&
			annot_regions->at(i_reg)->start >= profile_sites->at(i_site)->i_nuc)
		{
			// If the next site is after the start of the region, localization is successful.
			if(i_site < profile_sites->size() - 1 &&
				annot_regions->at(i_reg)->start < profile_sites->at(i_site+1)->i_nuc)
			{
				localized_site++;
				break;
			}
			i_site++;
		} // i_site loop.

		if(i_site < profile_sites->size())
		{
			// Dump the sites until the sites go out of range of the peak.
			char cur_peak_prof_fp[1000];
			sprintf(cur_peak_prof_fp, "%s/peak_%d_%d.sgr", dump_dir, annot_regions->at(i_reg)->start, annot_regions->at(i_reg)->end);
			FILE* f_cur_peak_prof = open_f(cur_peak_prof_fp, "w");

			// If file is not open, return.
			if(f_cur_peak_prof == NULL)
			{
				printf("Could not open %s for writing peak profile data.\n", cur_peak_prof_fp);
				return;
			}

			while(i_site < profile_sites->size() &&
				annot_regions->at(i_reg)->end > profile_sites->at(i_site)->i_nuc)
			{
				fprintf(f_cur_peak_prof, "%s\t%d\t%lf\n", annot_regions->at(i_reg)->chrom, profile_sites->at(i_site)->i_nuc, profile_sites->at(i_site)->height);

				i_site++;
			} // loop for processing current peak.

			fclose(f_cur_peak_prof);
		}
	} // i_ent loop.
}


t_signal_profile* load_sgr(char* sgr_fp)
{
	t_signal_profile* signal_profile = allocate_signal_profile();

	printf("Loading sgr file %s\n", sgr_fp);
	FILE* f_sgr = fopen(sgr_fp, "r");
	if(f_sgr == NULL)
	{
		printf("Could not open sgr file %s\n", sgr_fp);
		exit(0);
	}

	//vector<t_profile_site*>* prof_sites = new vector<t_profile_site*>();

	int cur_pos;
	double cur_height;
	char cur_chrom[100];
	while(fscanf(f_sgr, "%s %d %lf", cur_chrom, &cur_pos, &cur_height) == 3)
	{
		int i_chr = get_i_chr(signal_profile->chr_ids, cur_chrom);
		if(i_chr == signal_profile->chr_ids->size())
		{
			char* new_chr = new char[strlen(cur_chrom) + 2];
			strcpy(new_chr, cur_chrom);
			signal_profile->chr_ids->push_back(new_chr);
			i_chr = get_i_chr(signal_profile->chr_ids, cur_chrom);
			if(i_chr == signal_profile->chr_ids->size())
			{
				printf("Failed to correctly add the chromosome %s @ %s(%d)\n", cur_chrom, __FILE__, __LINE__);
				exit(0);
			}

			// Add a new profile site list for the new chromosome.
			signal_profile->profiles_per_chr->push_back(new vector<t_profile_site*>());
		}

		t_profile_site* new_prof_site = new t_profile_site();
		new_prof_site->i_nuc = cur_pos;
		new_prof_site->height = cur_height;

		// Add the new profile site to the list that is belongs to.
		signal_profile->profiles_per_chr->at(i_chr)->push_back(new_prof_site);
	}

	fclose(f_sgr);

	return(signal_profile);
}


/*
Following code loads whole bedGraph file, for all the chromosomes. bedGraph information is sortable, therefore there is no
problem loading the whole data at once.
*/
/*void load_bedGraph(char* bgr_fp, 
	t_signal_profile_list* profile_sites_per_chrom, 
	vector<char*>* chroms)*/
t_signal_profile* load_bedGraph(char* bgr_fp)
{
	t_signal_profile* signal_profile = allocate_signal_profile();

	printf("Loading bedGraph file %s\n", bgr_fp);
	FILE* f_bgr = fopen(bgr_fp, "r");
	if(f_bgr == NULL)
	{
		printf("Could not open bedGraph file %s\n", bgr_fp);
		exit(0);
	}

	// Get the chromosomes, this is very slow on large files.
	vector<vector<t_bedGraph_entry*>*>* bgr_entries_per_chrom = new vector<vector<t_bedGraph_entry*>*>();

	// Load the bedGraph entries, sort them.
	char cur_line[1000];
	int i_line = 0;
	while(x_fgets(cur_line, 1000, f_bgr))
	{
		//if((i_line % (100 * 1000)) == 0)
		//{
		//	fprintf(stderr, "%d             \r", i_line);
		//	i_line++;
		//}

		char cur_chrom[100];
		int start;
		int end;
		double signal;
		if(sscanf(cur_line, "%s %d %d %lf", cur_chrom, &start, &end, &signal) != 4)
		{
			printf("Could not parse the bedGraph line %s\n", cur_line);
			exit(0);
		}

		// Chromosome id check: Get the chromosome id.
		int i_chr = get_i_chr(signal_profile->chr_ids, cur_chrom);

		if(i_chr == signal_profile->chr_ids->size())
		{
			// Add this chromosome name to chromosome names and add a new list to the bgr_entries.
			char* new_chr_id = new char[strlen(cur_chrom) + 2];
			strcpy(new_chr_id, cur_chrom);
			signal_profile->chr_ids->push_back(new_chr_id);

			// Add the corresponding bgr entry list.
			bgr_entries_per_chrom->push_back(new vector<t_bedGraph_entry*>());

			printf("Added chromosome %s\n", cur_chrom);

			// Get the chromosome id.
			i_chr = get_i_chr(signal_profile->chr_ids, cur_chrom);

			// Sanity check.
			if(i_chr == signal_profile->chr_ids->size())
			{
				printf("Newly added chromosome %s was not found @ %s(%d)\n", cur_chrom, __FILE__, __LINE__);
				exit(0);
			}
		} // i_chr check.

		// Add this entry.
		t_bedGraph_entry* new_entry = new t_bedGraph_entry();
		new_entry->start = start;
		new_entry->end = end;
		new_entry->value = signal;
		bgr_entries_per_chrom->at(i_chr)->push_back(new_entry);
	} // file reading loop.

	// Allocate the profile site and bedGraph entry lists.
	for(int i_chr = 0; i_chr < signal_profile->chr_ids->size(); i_chr++)
	{
		signal_profile->profiles_per_chr->push_back(new vector<t_profile_site*>());
		printf("Loaded %d bgr entries at %s\n", bgr_entries_per_chrom->at(i_chr)->size(), signal_profile->chr_ids->at(i_chr));
	} // i_chr loop.

	// Sort the bedGraph entries.
	for(int i_chr = 0; i_chr < bgr_entries_per_chrom->size(); i_chr++)
	{
		printf("Generating profile sites for %s\n", signal_profile->chr_ids->at(i_chr));

		// Get the entries for current chromosome.
		vector<t_bedGraph_entry*>* bgr_entries = bgr_entries_per_chrom->at(i_chr);
		vector<t_profile_site*>* prof_sites = signal_profile->profiles_per_chr->at(i_chr);

		// Sort the bgr entries.
		sort(bgr_entries->begin(), bgr_entries->end(), sort_bedGraph_entries);

		 // Convert into profile sites: Add 0's between the necessary regions.
		int cur_index = 0;
		for(int i_ent = 0; i_ent < bgr_entries->size(); i_ent++)
		{
			// If the next region is farther away from from cur_index, assign zero to cur_index.
			// Note that the first 0 valued region, from 0 to beginning of the first region with value with 0 height, is skipped. 
			// This should also fix some problems in correlation computation.
			if(bgr_entries->at(i_ent)->start > cur_index &&
				cur_index > 0)
			{
				// Add a 0 height region between the cur_index and start of this bedGraph site.
				t_profile_site* new_site = new t_profile_site();

				// 0 valued region starts at the cur_index.
				new_site->i_nuc = cur_index;
				new_site->height = 0;
				//printf("Adding a 0 height region: %d\n", cur_index);
				prof_sites->push_back(new_site);
			}

			// Add the new profile_site starting at start.
			t_profile_site* new_site = new t_profile_site();
			new_site->i_nuc = bgr_entries->at(i_ent)->start;
			new_site->height = bgr_entries->at(i_ent)->value;

			prof_sites->push_back(new_site);

			// This is the next index without signal value, the next height assignment must start from this index.
			cur_index = bgr_entries->at(i_ent)->end;
		} // i_ent loop.

		// Add the last profile site.
		t_profile_site* new_site = new t_profile_site();

		// 0 valued region starts at the cur_index.
		new_site->i_nuc = cur_index;
		new_site->height = 0;
		prof_sites->push_back(new_site);

		// Delete the memory for bgr_entries in the current chromosome. 
		for(int i_bgr = 0; i_bgr < bgr_entries->size(); i_bgr++)
		{
			delete bgr_entries->at(i_bgr);
		} // i_bgr loop.

		delete bgr_entries;
	} // i_chr loop.

	delete(bgr_entries_per_chrom);

	printf("Loaded!\n");
	return(signal_profile);
}

/*
Determine all the peaks for this enrichment profile. 
*/
vector<t_peak_region*>** get_peaks_per_thresholds(vector<t_profile_site*>* profile_sites,
													int min_thresh,
													int max_thresh,
													int min_gap_bw_peaks,
													double scaling_factor)
{
	// Allocate the peak region lists, a list per threshold is allocated.	
	vector<t_peak_region*>** peak_regions_per_threshold = new vector<t_peak_region*>*[max_thresh - min_thresh + 1];
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		peak_regions_per_threshold[thresh - min_thresh] = new vector<t_peak_region*>();	
	} // thresh loop.	

    signed int* next_possible_peak_starts = new signed int[max_thresh - min_thresh + 1];
    bool* above_threshold = new bool[max_thresh - min_thresh + 1];
	bool* height_at_max = new bool[max_thresh - min_thresh + 1];
    for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
    {
            next_possible_peak_starts[thresh - min_thresh] = -10000;
            above_threshold[thresh - min_thresh] = false;
		height_at_max[thresh - min_thresh] = false;
    } // thresh loop.

    // Go over the profile and count the peaks.
    double height = 0;
	int cur_pos = 0;
	int n_sites = profile_sites->size();
	for(int i_site = 0; i_site < n_sites; i_site++)
    {
		// Update the current height: Note that there can be multiple ends at this position.
		height = scaling_factor * profile_sites->at(i_site)->height;
		cur_pos = profile_sites->at(i_site)->i_nuc;

        if(height < 0)
        {
                //fprintf(stderr, "Height below: %lf\n", height);
			    fprintf(stderr, "Height below 0: %.10f (%s, %d)\n", height, __FILE__, __LINE__);
                exit(0);
        }

		// At this point, the height is changed.
		// For each threhold that is lower than this height, this is either a new peak (if signal was below the threshold) otherwise
		// must check if the closest peak was at least 200 bps away.
		for(int thresh = min_thresh; 
			thresh <= max_thresh; 
			thresh++)
		{
			if(height >= (double)thresh)
			{
				if(above_threshold[thresh - min_thresh])
				{
					// The signal was above the threshold in the last signal update, thus it has to be above this threshold. Do not need to do anything at this point.

					// Process the maximum height for this peak, if it is necessary.
					if((int)height > peak_regions_per_threshold[thresh - min_thresh]->back()->max_height)
					{
						peak_regions_per_threshold[thresh - min_thresh]->back()->max_height = (int)height;
						peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_start_i = cur_pos;
						height_at_max[thresh - min_thresh] = true;
					}

					// Check if the height fell below the maximum.
					if(height_at_max[thresh - min_thresh] && (height < peak_regions_per_threshold[thresh - min_thresh]->back()->max_height))
					{
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_end_i = cur_pos - 1;
                                                height_at_max[thresh - min_thresh] = false;			
					}
				}
				else // Was below the threshold.
				{
					// check if this peak is far enough from the last peak.
					if(cur_pos > next_possible_peak_starts[thresh - min_thresh])
					{
						//printf("Peak @ %d\n", cur_pos);

						// Set the ending of the previous peak.
						if(peak_regions_per_threshold[thresh - min_thresh]->size() > 0)
						{
							peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
							// Also if the peak was in the max height, set it as the end of the peak.
							if(height_at_max[thresh - min_thresh])
							{
								peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_end_i = cur_pos - 1;
								height_at_max[thresh - min_thresh] = false;
							}
						}

						t_peak_region* new_peak_region = new t_peak_region(cur_pos, 0); // The end is not initialized, yet.
						peak_regions_per_threshold[thresh - min_thresh]->push_back(new_peak_region); // Add this peak region to the list of peaks.
						// Set the maximum height and position of the height for this peak when the peak is initiated.
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height = (int)height;
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_start_i = cur_pos;
                                                height_at_max[thresh - min_thresh] = true;
					}
				}

				// Set the new signal state as above threshold for this peak.
				above_threshold[thresh - min_thresh] = true;
			}
			else
			{
				// If the signal was above this threshold previously, then this is the end of a peak, mark the next possible peak start for this trheshold. Note that the maximum for the peak is not processed as the maximum height for a peak is always guaranteed to be at a position where the signal is higher than the threshold.
				if(above_threshold[thresh - min_thresh])
				{
					next_possible_peak_starts[thresh - min_thresh] = cur_pos + min_gap_bw_peaks-1;
				}
				above_threshold[thresh - min_thresh] = false;
			}
		} // thresh loop
	} // i_end loop.


	// Must check whether there is any peak that is hanging at the end. This is becaues a peak added when the next peak is starting.
	// Note that if there were no peaks in the whole chromosome, this loop should not be executed.
	for(int thresh = min_thresh; 
		thresh <= max_thresh;
		thresh++)
	{
		if(peak_regions_per_threshold[thresh - min_thresh]->size() > 0)
		{
			// If the last peak was below the threshold, this means that an end for this peak was found but the algorithm was looking for another peak to merge it with. Set the end for this peak.
			if(!above_threshold[thresh - min_thresh])
			{
				peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
			}
			else
			{
				 peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = cur_pos;				
			}
		}
	} // thresh loop

	return(peak_regions_per_threshold);	
}

//void load_sgr(char* sgr_fp, 
//				t_signal_profile_list* profile_sites_per_chroms,
//				vector<char*>* chr_ids)
//{
//	printf("Loading sgr file %s\n", sgr_fp);
//	FILE* f_sgr = fopen(sgr_fp, "r");
//	if(f_sgr == NULL)
//	{
//		printf("Could not open sgr file %s\n", sgr_fp);
//		exit(0);
//	}
//
//	//vector<t_profile_site*>* prof_sites = new vector<t_profile_site*>();
//
//	int cur_pos;
//	double cur_height;
//	char cur_chrom[100];
//	while(fscanf(f_sgr, "%s %d %lf", cur_chrom, &cur_pos, &cur_height) == 3)
//	{
//		int i_chr = get_i_chr(chr_ids, cur_chrom);
//		if(i_chr == chr_ids->size())
//		{
//			char* new_chr = new char[strlen(cur_chrom) + 2];
//			strcpy(new_chr, cur_chrom);
//			chr_ids->push_back(new_chr);
//			i_chr = get_i_chr(chr_ids, cur_chrom);
//			if(i_chr == chr_ids->size())
//			{
//				printf("Failed to correctly add the chromosome %s @ %s(%d)\n", cur_chrom, __FILE__, __LINE__);
//				exit(0);
//			}
//
//			// Add a new profile site list for the new chromosome.
//			profile_sites_per_chroms->push_back(new vector<t_profile_site*>());
//		}
//
//		t_profile_site* new_prof_site = new t_profile_site();
//		new_prof_site->i_nuc = cur_pos;
//		new_prof_site->height = cur_height;
//
//		// Add the new profile site to the list that is belongs to.
//		profile_sites_per_chroms->at(i_chr)->push_back(new_prof_site);
//	}
//
//	fclose(f_sgr);
//}

vector<t_profile_site*>* get_delta_profile_per_profile(vector<t_profile_site*>* profile)
{
	sort(profile->begin(), profile->end(), sort_profile_sites);

	// Sort the profile sites.
	vector<t_profile_site*>* delta_profile = new vector<t_profile_site*>();

	double prev_height = 0.0;
	int prev_pos = 0;

	for(int i_site = 0; i_site < profile->size(); i_site++)
	{
		//if(profile->at(i_site)->height != prev_height)
		//{
			// Add the fragment ends at this position: If the new height is larger, a start is required.
			double cur_dh = (profile->at(i_site)->height - prev_height);

			t_profile_site* new_site = new t_profile_site();
			new_site->i_nuc = profile->at(i_site)->i_nuc;
			new_site->height = cur_dh;

			delta_profile->push_back(new_site);

			//for(int h = prev_height; h != profile->at(i_site)->height; h += cur_dh)
			//{
				//t_frag_end* new_end = new t_frag_end();
				//new_end->pos = profile->at(i_site)->i_nuc;
				//new_end->side = end_type;

			//	frag_ends->push_back(new_end);
			//} // h loop
		//}
		//else
		//{
		//	printf("Skipping %d-%lf (%lf)\n", profile->at(i_site)->i_nuc, profile->at(i_site)->height, prev_height);
		//}

		// Update previous height and position.
		prev_height = profile->at(i_site)->height;
		prev_pos = profile->at(i_site)->i_nuc;
	} // i_site loop.

	return(delta_profile);
}

vector<t_profile_site*>* get_profile_per_delta_profile(vector<t_profile_site*>* delta_profile)
{
	// Sort delta profile.
	sort(delta_profile->begin(), delta_profile->end(), sort_profile_sites);

	// Allocate profile.
	vector<t_profile_site*>* profile_sites = new vector<t_profile_site*>();

	double height = 0;
	int cur_pos = 1;
    for(int i_end = 0; i_end < delta_profile->size(); i_end++)
    {
        // Update the current height: Note that there can be multiple ends at this position.
		height += delta_profile->at(i_end)->height;
        cur_pos = delta_profile->at(i_end)->i_nuc;

        // Now add all the other ends at this position.
        while((i_end + 1) < delta_profile->size() &&
				cur_pos == delta_profile->at(i_end+1)->i_nuc)
        {
                i_end++;
                height += delta_profile->at(i_end)->height;
        }

		// The newly added region overlaps with the requested region, dump the previous region and height, since 
		// the region corresponds to the height of the revious pos.
        //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);
		t_profile_site* new_site = new t_profile_site();
		new_site->i_nuc = cur_pos;
		new_site->height = height;
		profile_sites->push_back(new_site);
    } // i_frag

	return(profile_sites);
}

void delete_signal_profile(t_signal_profile* signal_profile)
{
	for(int i_chr = 0; i_chr < signal_profile->profiles_per_chr->size(); i_chr++)
	{
		delete_profile_sites(signal_profile->profiles_per_chr->at(i_chr));
		delete [] signal_profile->chr_ids->at(i_chr);
	} // i_chr loop.

	delete(signal_profile->profiles_per_chr);
	delete(signal_profile->chr_ids);
	delete(signal_profile);
}

bool sort_profile_sites(t_profile_site* entry1, t_profile_site* entry2)
{
	return(entry1->i_nuc < entry2->i_nuc);
}

vector<t_annot_region*>* load_BED_with_avg_features(char* local_means_fp, int block_length)
{
	vector<t_annot_region*>* local_mean_info = load_BED_with_line_information(local_means_fp);

	//vector<vector<double>*>* local_means_per_chr = new vector<vector<double>*>();
	vector<t_annot_region*>* local_mean_regions = new vector<t_annot_region*>();

	for(int i_reg = 0; i_reg < local_mean_info->size(); i_reg++)
	{
		// Initialize the region information.
		t_annot_region* new_chr_region = new t_annot_region();
		new_chr_region->chrom = new char[strlen(local_mean_info->at(i_reg)->chrom) + 2];
		strcpy(new_chr_region->chrom, local_mean_info->at(i_reg)->chrom);
		new_chr_region->start = local_mean_info->at(i_reg)->start;
		new_chr_region->end = local_mean_info->at(i_reg)->end;
		new_chr_region->strand = local_mean_info->at(i_reg)->strand;

		char* cur_line = (char*)(local_mean_info->at(i_reg)->data);

		t_string* cur_line_str = new t_string(cur_line);

		// Tokenize the current line, read the local averages: First 6 entries are coordinates, chromsomes and etc.
		t_string_tokens* line_str_tokens = cur_line_str->tokenize_by_chars(" \t");

		vector<double>* cur_chr_means = new vector<double>();
		for(int i_tok = 6; i_tok < line_str_tokens->size(); i_tok++)
		{
			double current_local_mean = atof(line_str_tokens->at(i_tok)->str());
			current_local_mean /= block_length;
			cur_chr_means->push_back(current_local_mean);
		} // i_tok loop.

		new_chr_region->data = (void*)cur_chr_means;

		// Add the new region.
		local_mean_regions->push_back(new_chr_region);

		t_string::clean_tokens(line_str_tokens);
		delete(cur_line_str);
	} // i_reg loop.

	// Free memory.
	delete_annot_regions_with_line_information(local_mean_info);

	return(local_mean_regions);
}









