#include <stdio.h>
#include <stdlib.h>
#include "../../../lib/chromosome/chromosome.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "../../../lib/database/mapped_read/fragment.h"
#include <algorithm>
#include "simulated_background.h"
#include <math.h>
#include <string.h>

double get_input_scaling_factor(int n_meg_wins,
				vector<t_fragment*>* chip_seq_fragments,
				vector<t_fragment*>* control_fragments,
				double P_f, 
				int enrichment_fragment_length)
{
	printf("Computing input scaling for each chromosome.\n");

	int chip_seq_background_cnt = 0;
	int control_cnt = 0;

/*
	// Load next 250 lines in mappability map file.
	for(int i_win = 0; i_win < 250; i_win++)
	{
		int read_i_win = 0;
		fscanf(f_mapp_map, "%*s %d %d\n", &read_i_win, &uniquely_mappable_base_cnt_per_win[i_win]);
		if(read_i_win != i_win)
		{
			printf("Inconsistent mappability file (%d, %d) @ %s(%d)\n", read_i_win, i_win, __FILE__, __LINE__);
			exit(0);
		}
	} // i_win loop.
	fprintf(stderr, " [done]\n");
*/

	vector<int>* chip_seq_frag_cnt = new vector<int>();
	vector<int>* control_frag_cnt = new vector<int>();

	// For each of the chromosomes, load the fragment counts.
	//for(int i_chr = 21; i_chr < chr_fps->size(); i_chr++)
/*
	for(int i = 0; i < 20; i++)
	{
		printf("%d. fragment: %d\n", i, cur_chr_fragments->at(i)->base_index);
	}

	fprintf(stderr, "Building enrichment profile: ");
	// Generate the enrichment profile for this chromosome.
	t_enrichment_profile* cur_chr_enrichment_profile = new t_enrichment_profile(cur_chr->l_chromosome(), cur_chr_chip_seq_fragments, enrichment_fragment_length);
	fprintf(stderr, "[done]\n");

	for(int i = 0; i < 3000; i++)
	{
		if(cur_chr_enrichment_profile->enrichment_counts[i] > 0)
		{
			printf("%d: %d\n", i, cur_chr_enrichment_profile->enrichment_counts[i]);
		}
	} // i loop.

        FILE* f_chip_seq_frags = open_f("chip_seq_frags.txt", "w");
        for(int i = 0; i < cur_chr_chip_seq_fragments->size(); i++)
        {
                printf("%d: %d\n", i, sample_starts->ps[i]);
                fprintf(f_chip_seq_frags, "%d: %d\n", i, cur_chr_chip_seq_fragments->at(i)->base_index);
        }
        fclose(f_chip_seq_frags);

        FILE* f_control_frags = open_f("control_frags.txt", "w");
        printf("Control fragments:\n");
        for(int i = 0; i < cur_chr_control_fragments->size(); i++)
        {
                //printf("%d: %d\n", i, control_starts->ps[i]);
                fprintf(f_control_frags, "%d: %d\n", i, cur_chr_control_fragments->at(i)->base_index);
        }
        fclose(f_control_frags);
*/
/*
	for(int i = 0; i < 20; i++)
	{
		printf("%d. fragment: %d\n", i, cur_chr_control_fragments->at(i)->base_index);
	}
*/
	//getc(stdin);
	// Now, determine the threshold on enrichment counts  for which P_f is satisfied, this will be used in the next loop to filter out
	// the (scaling) windows 

	int win_size = 10 * K_BASE;

	// Count the numbers per window. 
	int n_scaling_win = n_meg_wins * MEG_BASE / win_size;

	//int chip_seq_background_cnt = 0;
	//int control_cnt = 0;

	// Exploit the fact that the fragments are sorted with respect to base indices in the lists:
	// Keep track of the current fragment index for both fragment lists.
	int cur_chr_fragment_list_base_i = 0;
	int cur_chr_control_fragment_list_base_i = 0;

	//FILE* f_dps = open_f("dps.txt", "w");

	int n_processed_wins = 0;

	fprintf(stderr, "Computing scaling factor over %d windows (%d, %d).\n", n_scaling_win, chip_seq_fragments->size(), control_fragments->size());
	for(int i_scaling_win = 0; i_scaling_win < n_scaling_win; i_scaling_win++)
	{
		if(i_scaling_win % 100 == 0)
		{
			//fprintf(stderr, "%d(%d)            \r", i_scaling_win, n_scaling_win);
		}

		int min_i_nuc = win_size * i_scaling_win;
		int max_i_nuc = (i_scaling_win + 1) * win_size;

		// Firstly, ensure that this scaling window does not contain any significant nucleotide peaks in it.
		// If so, count the reads in this region from the fragment list.
		// The check is done by checking if there is any peak in the current window whose count is 
		// higher than the threshold computed above to satisfy P_f to exclude certain percentage of nucleotide peaks.
		// Note that the fragments are sorted with respect to their positions on the chromosome.
		// This is done by determining if the enrichment profile of ChIP-Seq data has a peak whose enrichment count
		// is higher than the threshold determined using P_f.


		// Now, if the window passed the test, process it: Count all the reads in the window.
		int cur_win_chip_seq_bgrnd_cnt = 0;
		int cur_win_control_cnt = 0;
		if(true)
		{
			// [PeakSeq compatilibity]
			int last_added_chip_i = -1;
			//printf("At %d (%d-%d)\n", cur_chr_data->chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index, min_i_nuc, max_i_nuc);

			// Find the index that is in the current window limits.
            while(cur_chr_fragment_list_base_i < chip_seq_fragments->size() &&
                    chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index < min_i_nuc)
            {
				cur_chr_fragment_list_base_i++;
			}

			while(cur_chr_fragment_list_base_i < chip_seq_fragments->size() &&
				chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index >= min_i_nuc &&
				chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index < max_i_nuc)
            {
				//chip_seq_background_cnt++;
			
				//printf("At %d (%d)\n", cur_chr_data->chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index, last_added_chip_i);		
				if(last_added_chip_i != chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index)
				{
					last_added_chip_i = chip_seq_fragments->at(cur_chr_fragment_list_base_i)->base_index;

					cur_win_chip_seq_bgrnd_cnt++;
				}

				cur_chr_fragment_list_base_i++;
            } // i_nuc loop.

			// [PeakSeq compatilibity]	

			// Note that followin loop is necesasary but taken out to comply with
            // PeakSeq.

/*
			// Find the index that is in the current window limits.
            while(cur_chr_control_fragment_list_base_i < cur_chr_data->control_fragments->size() &&
                    cur_chr_data->control_fragments->at(cur_chr_control_fragment_list_base_i)->base_index < min_i_nuc)
            {
				cur_chr_control_fragment_list_base_i++;
			}
*/
		int last_added_control_i = -1;

        while(cur_chr_control_fragment_list_base_i < control_fragments->size() &&
                //cur_chr_data->control_fragments->at(cur_chr_control_fragment_list_base_i)->base_index >= min_i_nuc &&
                control_fragments->at(cur_chr_control_fragment_list_base_i)->base_index < max_i_nuc)
        {
                //control_cnt++;

				if(last_added_control_i != control_fragments->at(cur_chr_control_fragment_list_base_i)->base_index)
				{
					cur_win_control_cnt++;	
					last_added_control_i = control_fragments->at(cur_chr_control_fragment_list_base_i)->base_index;
				}

				cur_chr_control_fragment_list_base_i++;
            } // i_nuc loop.

			//if(cur_win_chip_seq_bgrnd_cnt > 0 && cur_win_control_cnt > 0)
			if(cur_win_chip_seq_bgrnd_cnt > 0 && cur_win_control_cnt > 0)
			{
               	//fprintf(stderr, "%d %d\n", cur_win_chip_seq_bgrnd_cnt, cur_win_control_cnt);
				chip_seq_frag_cnt->push_back(cur_win_chip_seq_bgrnd_cnt);
				control_frag_cnt->push_back(cur_win_control_cnt);

                //FILE* f_bin_cnts = open_f("bin_cnts.txt", "a");
                //fprintf(f_bin_cnts, "bin %d, control count: %d, sample count: %d\n", i_scaling_win, cur_win_control_cnt, cur_win_chip_seq_bgrnd_cnt);
                //fclose(f_bin_cnts);

				n_processed_wins++;
			}
		}
	} // i_scaling_win loop.

        // Compute the slope.
	double slope = PeakSeq_slope(control_frag_cnt, chip_seq_frag_cnt, n_processed_wins);
	
        double r = pearson_correlation(control_frag_cnt, chip_seq_frag_cnt);

        printf("Slope: %lf\nPearson Coeff: %lf\n", slope, r);

	//return(scaling_factor_per_chromosome);
	return(slope);
}

double PeakSeq_slope(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt, int n_processed_wins)
{
        // Compute the slope from 2 DOF (slope) LLS regression.
        double num = 0.0f;
        double den = 0.0f;

        double sx = 0.0f;
        double sy = 0.0f;
        double sxx = 0.0f;
        double sxy = 0.0f;

        for(int i = 0; i < control_frag_cnt->size(); i++)
        {
                sx += control_frag_cnt->at(i);
                sy += chip_seq_frag_cnt->at(i);
                sxx += control_frag_cnt->at(i) * control_frag_cnt->at(i);
                sxy += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
        } // i loop

	//int n_processed_wins = control_frag_cnt->size();

	printf("n_processed_wins = %d (%lf, %lf, %lf, %lf)\n", n_processed_wins, sx, sy, sxx, sxy);

        num = n_processed_wins * sxy - sx * sy;
        den = n_processed_wins * sxx - sx * sx;

	return(num / den);
}

double slope(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt)
{
        // Compute the slope from 1 DOF (slope) LLS regression.
        double num = 0.0f;
        double den = 0.0f;
        for(int i = 0; i < control_frag_cnt->size(); i++)
        {
                num += control_frag_cnt->at(i) * chip_seq_frag_cnt->at(i);
                den += control_frag_cnt->at(i) * control_frag_cnt->at(i);
        } // i loop

	return(num / den);
}


// Copied from http://www.vias.org/tmdatanaleng/cc_corr_coeff.html
// Also check http://www.statsdirect.com/help/regression_and_correlation/sreg.htm
double pearson_correlation(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt)
{
	double num = 0.0f;
	double den = 0.0f;

	double control_mean = 0.0f;
	for(int i_frag = 0; i_frag < control_frag_cnt->size(); i_frag++)
	{
		control_mean += control_frag_cnt->at(i_frag);
	} // i_frag loop.
	control_mean /= control_frag_cnt->size();

        double chip_seq_mean = 0.0f;
        for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
        {
                chip_seq_mean += chip_seq_frag_cnt->at(i_frag);
        } // i_frag loop.
	chip_seq_mean /= chip_seq_frag_cnt->size();

        for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
        {
		num += (control_frag_cnt->at(i_frag) - control_mean) * (chip_seq_frag_cnt->at(i_frag) - chip_seq_mean);
	} // i_frag loop.
	
	double control_var = 0.0f;

        for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
        {
		control_var += (control_frag_cnt->at(i_frag) - control_mean) * (control_frag_cnt->at(i_frag) - control_mean);
	}	
	control_var = pow(control_var, .5);

	double chip_seq_var = 0.0f;
        for(int i_frag = 0; i_frag < chip_seq_frag_cnt->size(); i_frag++)
        {
		chip_seq_var += (chip_seq_frag_cnt->at(i_frag) - chip_seq_mean) * (chip_seq_frag_cnt->at(i_frag) - chip_seq_mean);
	}
	chip_seq_var = pow(chip_seq_var, .5);

	double r = num / (control_var * chip_seq_var);

	return(r);
}



