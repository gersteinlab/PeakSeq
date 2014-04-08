#include <stdio.h>
#include <stdlib.h>
#include "chip_seq_chr_data.h"
#include "../../../lib/database/mapped_read/mapped_read_tools.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "simulated_background.h"
#include "../../../lib/chromosome/chromosome.h"
#include "compare_signal_tracks.h"
#include "../../../lib/database/mapped_read/fragment.h"
#include "peakseq.h"
#include "peak_region.h"
#include "../../../lib/utils/file/utils.h"
#include <string.h>

/*
class t_chip_seq_chr_data
{
public:
        // The id is used to form the file names associated with this chromosome.
        t_chip_seq_chr_data(char* chr_id, char* mappability_map_fp);
        ~t_chip_seq_chr_data();

        t_chromorome* chr; // Can be NULL.

        // The set of chip seq and control fragments in the
        vector<t_fragment*>* chip_seq_fragments;
        vector<t_fragment*>* control_fragment;

        // The enrichment profiles.
        t_enrichment_profile* chip_seq_enr_prof;
        t_enrichment_profile* control_enr_prof;

        // Number of uniquely mappable nucleotide counts per megabase window.
        vector<int*>* n_uniquely_mappable_nucs_per_meg_win;

        // Number of megabase windows in this chromosome: Use the
        int n_meg_wins();
};
*/

t_chip_seq_chr_data::t_chip_seq_chr_data(vector<char*>* chip_seq_reads_data_dirs, 
					vector<char*>* input_reads_data_dirs, 
					char* _chr_id, 
					char* mappability_map_fp, 
					int enrichment_fragment_length)
{
    FILE* f_mapp_map = open_f(mappability_map_fp, "r");
    if(f_mapp_map == NULL)
    {
            printf("Could not open mappability map.\n");
            exit(0);
    }

	// Copy chromosome id.
	this->chr_id = new char[strlen(_chr_id) + 2];
	strcpy(this->chr_id, _chr_id);

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
				printf("The ordering in Mapability map file must be fixed.\n");
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
		printf("Did not read any mappable nucleotide for chromosome %s. Ensure that the id's in chromosome id list are consistent with the ids in %s.\n", this->chr_id, mappability_map_fp);
		//exit(0);
		this->cur_chr_chip_seq_fore_frags = NULL;
		this->cur_chr_chip_seq_rev_frags = NULL;
		this->chip_seq_fragments = NULL;
		this->control_fragments = NULL;
		this->chip_seq_enr_prof = NULL;
		this->control_enr_prof = NULL;
	}
	
	// Load the fragments.
	fprintf(stderr, "Loading ChIP-Seq fragments from %d directories: ", chip_seq_reads_data_dirs->size());
//        vector<t_fragment*>* cur_chr_chip_seq_fore_frags = new vector<t_fragment*>();
//        vector<t_fragment*>* cur_chr_chip_seq_rev_frags = new vector<t_fragment*>();
	this->cur_chr_chip_seq_fore_frags = new vector<t_fragment*>();
	this->cur_chr_chip_seq_rev_frags = new vector<t_fragment*>();

    // Load all the data from each directory.
    for(int i_chip_seq_dir = 0; i_chip_seq_dir < chip_seq_reads_data_dirs->size(); i_chip_seq_dir++)
    {	
		char mapped_chip_seq_reads_fp[10000];
		sprintf(mapped_chip_seq_reads_fp, "%s/%s_mapped_reads.bin", chip_seq_reads_data_dirs->at(i_chip_seq_dir), this->chr_id); 
		//printf("Loading %s\n", mapped_chip_seq_reads_fp);
		load_fragments_binary(mapped_chip_seq_reads_fp, cur_chr_chip_seq_fore_frags, cur_chr_chip_seq_rev_frags);
	} // i_chip_seq_dir loop.

	this->chip_seq_fragments = forwardize_combine_sort_fore_rev_strand_frags(cur_chr_chip_seq_fore_frags, cur_chr_chip_seq_rev_frags, enrichment_fragment_length);
	//delete_fragments(cur_chr_chip_seq_fore_frags);
	//delete_fragments(cur_chr_chip_seq_rev_frags);

	printf("[done] Loaded %d fragments.\n", this->chip_seq_fragments->size());

	fprintf(stderr, "Building enrichment profile: ");
	// Generate the enrichment profile for this chromosome.
	this->chip_seq_enr_prof = new t_enrichment_profile( this->chr_id, this->chip_seq_fragments, enrichment_fragment_length);
	fprintf(stderr, "[done]\n");

/*
	char enr_prof_fp[1000];
	sprintf(enr_prof_fp, "%s_chip_seq_enr_prof.txt", this->chr_id);
	this->chip_seq_enr_prof->dump_profile(enr_prof_fp, 1, 100000000);
*/

	// Load control fragments and generate the profile
	fprintf(stderr, "Loading control fragments from %d directories: ", input_reads_data_dirs->size());
        vector<t_fragment*>* cur_chr_fore_control_frags = new vector<t_fragment*>();
        vector<t_fragment*>* cur_chr_rev_control_frags = new vector<t_fragment*>();

	for(int i_input_dir = 0; i_input_dir < input_reads_data_dirs->size(); i_input_dir++)
	{
		char mapped_control_reads_fp[10000];
		sprintf(mapped_control_reads_fp, "%s/%s_mapped_reads.bin", input_reads_data_dirs->at(i_input_dir), this->chr_id);
		//printf("Loading %s\n", mapped_control_reads_fp);
		//sprintf(mapped_control_reads_fp, "parsed_txt_reads/%s_control_mapped_reads.txt", chr_id);

	        load_fragments_binary(mapped_control_reads_fp, cur_chr_fore_control_frags, cur_chr_rev_control_frags);
	        //load_fragments(NULL, mapped_control_reads_fp, cur_chr_fore_control_frags, cur_chr_rev_control_frags);
	} // i_input_dir loop.

	this->control_fragments = forwardize_combine_sort_fore_rev_strand_frags(cur_chr_fore_control_frags, cur_chr_rev_control_frags, enrichment_fragment_length);

/*
    for(int i_frag = 0; i_frag < 100; i_frag++)
    {
            printf("%d\n", this->control_fragments->at(i_frag)->base_index);
    }

    exit(0);
*/

	delete_fragments(cur_chr_fore_control_frags);
	delete_fragments(cur_chr_rev_control_frags);

	printf("[done] Loaded %d fragments.\n", this->control_fragments->size());

	fprintf(stderr, "Building enrichment profile: ");
	// Generate the enrichment profile for this chromosome.
	this->control_enr_prof = new t_enrichment_profile(this->chr_id, this->control_fragments, enrichment_fragment_length);
	fprintf(stderr, "[done]\n");
	
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

int t_chip_seq_chr_data::n_meg_wins()
{
        return(this->n_uniquely_mappable_nucs_per_meg_win->size());
}
