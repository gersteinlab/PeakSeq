#include <stdio.h>
#include <stdlib.h>
#include "chip_seq_chr_data.h"
#include "mapped_read_tools.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "simulated_background.h"
#include "compare_signal_tracks.h"
#include "peakseq.h"
#include "peak_region.h"
#include "utils.h"
#include <string.h>
#include <math.h>

bool __DUMP_CHIP_SEQ_CHR_DATA_MSGS__ = false;

t_chip_seq_chr_data::t_chip_seq_chr_data(vector<char*>* chip_seq_reads_data_dirs, 
					vector<char*>* input_reads_data_dirs, 
					char* _chr_id, 
					char* mappability_map_fp, 
					int enrichment_mapped_fragment_length)
{
	// Copy chromosome id.
	this->chr_id = new char[strlen(_chr_id) + 2];
	strcpy(this->chr_id, _chr_id);

	// Load the fragments.
if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
	fprintf(stderr, "Loading ChIP-Seq fragments from %ld directories: ", chip_seq_reads_data_dirs->size());

	this->cur_chr_chip_seq_fore_frags = new vector<t_mapped_fragment*>();
	this->cur_chr_chip_seq_rev_frags = new vector<t_mapped_fragment*>();

    // Load all the data from each directory.
    for(int i_chip_seq_dir = 0; i_chip_seq_dir < (int)chip_seq_reads_data_dirs->size(); i_chip_seq_dir++)
    {	
		char mapped_chip_seq_reads_fp[10000];
		sprintf(mapped_chip_seq_reads_fp, "%s/%s_mapped_reads.txt", chip_seq_reads_data_dirs->at(i_chip_seq_dir), this->chr_id); 
		load_fragments(mapped_chip_seq_reads_fp, this->cur_chr_chip_seq_fore_frags, this->cur_chr_chip_seq_rev_frags);
	} // i_chip_seq_dir loop.

	this->chip_seq_fragments = forwardize_combine_sort_fore_rev_strand_frags(cur_chr_chip_seq_fore_frags, cur_chr_chip_seq_rev_frags, enrichment_mapped_fragment_length);
	delete_fragments(cur_chr_chip_seq_fore_frags);
	delete_fragments(cur_chr_chip_seq_rev_frags);

if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
{
	fprintf(stderr, "Loaded %ld fragments.\n", this->chip_seq_fragments->size());
	fprintf(stderr, "Building enrichment profile: ");
}
	// Generate the enrichment profile for this chromosome.
	this->chip_seq_enr_prof = new t_enrichment_profile(this->chr_id, this->chip_seq_fragments, enrichment_mapped_fragment_length);

	// Load control fragments and generate the profile.
if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
{
	fprintf(stderr, "Loading control fragments from %ld directories: ", input_reads_data_dirs->size());
}

    vector<t_mapped_fragment*>* cur_chr_fore_control_frags = new vector<t_mapped_fragment*>();
    vector<t_mapped_fragment*>* cur_chr_rev_control_frags = new vector<t_mapped_fragment*>();

	for(int i_input_dir = 0; i_input_dir < (int)input_reads_data_dirs->size(); i_input_dir++)
	{
		char mapped_control_reads_fp[10000];
		//sprintf(mapped_control_reads_fp, "%s/%s_mapped_reads.bin", input_reads_data_dirs->at(i_input_dir), this->chr_id);
		sprintf(mapped_control_reads_fp, "%s/%s_mapped_reads.txt", input_reads_data_dirs->at(i_input_dir), this->chr_id);
		//printf("Loading %s\n", mapped_control_reads_fp);
		//sprintf(mapped_control_reads_fp, "parsed_txt_reads/%s_control_mapped_reads.txt", chr_id);

	    //load_fragments_binary(mapped_control_reads_fp, cur_chr_fore_control_frags, cur_chr_rev_control_frags);
		load_fragments(mapped_control_reads_fp, cur_chr_fore_control_frags, cur_chr_rev_control_frags);
	    //load_fragments(NULL, mapped_control_reads_fp, cur_chr_fore_control_frags, cur_chr_rev_control_frags);
	} // i_input_dir loop.

	this->control_fragments = forwardize_combine_sort_fore_rev_strand_frags(cur_chr_fore_control_frags, cur_chr_rev_control_frags, enrichment_mapped_fragment_length);

	delete_fragments(cur_chr_fore_control_frags);
	delete_fragments(cur_chr_rev_control_frags);

if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
{
	printf("[done] Loaded %ld fragments.\n", this->control_fragments->size());
	fprintf(stderr, "Building enrichment profile: ");
}
	// Generate the enrichment profile for this chromosome.
	this->control_enr_prof = new t_enrichment_profile(this->chr_id, this->control_fragments, enrichment_mapped_fragment_length);

	this->load_mapability_map_file(mappability_map_fp);
}

void t_chip_seq_chr_data::load_mapability_map_file(char* mappability_map_fp)
{
    FILE* f_mapp_map = fopen(mappability_map_fp, "r");

	// Check the existence of the mapability map file.
    if(f_mapp_map == NULL)
    {
		// This case is tricky. Setup the mapability map as all the windows being totally mappable.
		int n_wins = (int)(floor((double)(this->n_nucleotides()) / MEG_BASE))+1;

if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
		fprintf(stderr, "Setting up the full mapability map file for %d nucleotides (%d windows.)", this->n_nucleotides(), n_wins);

		this->n_uniquely_mappable_nucs_per_meg_win = new vector<int>();

		// All the windows are totally mappable.
		for(int i_win = 0; i_win < n_wins; i_win++)
		{
			this->n_uniquely_mappable_nucs_per_meg_win->push_back(MEG_BASE);
		} // i_win loop.
    }
	else
	{
if(__DUMP_CHIP_SEQ_CHR_DATA_MSGS__)
		fprintf(stderr, "Loading uniquely mappable nucleotide counts from %s: ", mappability_map_fp);

		this->n_uniquely_mappable_nucs_per_meg_win = new vector<int>();

		char cur_chr_id[100];
		int read_i_win = 0;

		int cur_n_unique_base_cnt = 0;

		int last_read_i_win = -1;
		while(fscanf(f_mapp_map, "%s %d %d\n", cur_chr_id, &read_i_win, &cur_n_unique_base_cnt) == 3)
		{
			// Save the read number if the id matches with the current chromosome id.
			if(strcmp(cur_chr_id, this->chr_id) == 0)
    		{
					// If there is a discontinuity in the window indices, fill up the skipped windows with 0 uniquely mappable nucleotides.
					if(read_i_win != last_read_i_win + 1)
					{
						if(read_i_win < last_read_i_win)
						{
							printf("The ordering in Mapability map file must be fixed; %d proceeds %d\n", read_i_win, last_read_i_win);
							exit(0);
						}

						for(int i_win = last_read_i_win+1; i_win <= read_i_win-1; i_win++)
						{
							this->n_uniquely_mappable_nucs_per_meg_win->push_back(0);
						} // i_win loop.
					}

					//printf("%s (%d): %d mappable nucs.\n", chr_id, read_i_win, cur_n_unique_base_cnt);
					this->n_uniquely_mappable_nucs_per_meg_win->push_back(cur_n_unique_base_cnt);

    				// This is here to make sure that windows follow up correctly.
					last_read_i_win = read_i_win;
			}
		} // i_win loop.
		//fprintf(stderr, " %d windows [done]\n", n_wins);

		// Sanity check, check if the mappability map is read. This is important to ensure that the id's in Mapability map are consistent with the 
		// ids from the id list file.
		if( this->n_uniquely_mappable_nucs_per_meg_win->size() == 0)
		{
			fprintf(stderr, "Did not read any mappable nucleotide for chromosome %s. Ensure that the id's in chromosome id list are consistent with the ids in %s.\n", this->chr_id, mappability_map_fp);
			//exit(0);
			this->cur_chr_chip_seq_fore_frags = NULL;
			this->cur_chr_chip_seq_rev_frags = NULL;
			this->chip_seq_fragments = NULL;
			this->control_fragments = NULL;
			this->chip_seq_enr_prof = NULL;
			this->control_enr_prof = NULL;
		}

		fclose(f_mapp_map);
	} // Check for existence of mapability map file.
}

t_chip_seq_chr_data::~t_chip_seq_chr_data()
{
        // The set of chip seq and control fragments in the
		if(this->chip_seq_fragments != NULL)
			delete_fragments(this->chip_seq_fragments);

		if(this->control_fragments != NULL)
			delete_fragments(this->control_fragments);

        // The enrichment profiles.
		if(this->chip_seq_enr_prof != NULL)
			delete(this->chip_seq_enr_prof);

		if(this->control_enr_prof != NULL)
			delete(this->control_enr_prof);

	delete [] this->chr_id;

	delete(this->n_uniquely_mappable_nucs_per_meg_win);
}

int t_chip_seq_chr_data::n_nucleotides()
{
	// Use the fragment list to determine the number of nucleotides.
	return(this->chip_seq_fragments->back()->base_index + this->chip_seq_fragments->back()->sequenced_fragment_length);
}

int t_chip_seq_chr_data::n_meg_wins()
{
        return(this->n_uniquely_mappable_nucs_per_meg_win->size());
}


