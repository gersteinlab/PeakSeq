#include <stdio.h>
#include <stdlib.h>
#include "simulated_background.h"
#include "enrichment_profile.h"
#include "mapped_read_tools.h"
#include "rng.h"
#include "seed_manager.h"
#include <math.h>
#include "chip_seq_chr_data.h"
#include "xlog_math.h"
#include <algorithm>

bool __DUMP_POISSON_BCKGRND_MSGS__ = false;

//int n_rand_gens = 0;
//FILE* f_rand = NULL;

double* get_poisson_thresholds(t_chip_seq_chr_data* cur_chr_data,
								double target_fdr,
								int enrichment_mapped_fragment_length,
								int min_thresh,
								int max_thresh)
{
	if(__DUMP_POISSON_BCKGRND_MSGS__)
	{
		fprintf(stderr, "Computing Poisson based thresholds %d,..., %d\n", min_thresh, max_thresh);
	}

	int n_wins = cur_chr_data->n_meg_wins();

	double* thresholds_per_win = new double[n_wins+2];
	for(int i_win = 0; i_win <= n_wins; i_win++)
	{
		//thresholds_per_win[i_win] = max_thresh; // Ensure that no peaks are selected at the initialization, this automatically handles the windows for which there is no data.
		thresholds_per_win[i_win] = 0;
	}

	// Count the number of reads per window.
	int* n_reads_per_window = get_n_reads_per_window(cur_chr_data->n_meg_wins(), cur_chr_data->chip_seq_fragments);

	// For each window, go over the threholds to find the threholds for which the probability is above 95%.
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		if(n_reads_per_window[i_win] == 0)
		{
			if(__DUMP_POISSON_BCKGRND_MSGS__)
				printf("There are no mapped reads in window %d\n", i_win);
			//getc(stdin);
		}

		if(cur_chr_data->n_uniquely_mappable_nucs_per_meg_win->at(i_win) == 0)
		{
			if(__DUMP_POISSON_BCKGRND_MSGS__)
				printf("There are no uniquely mappable nucleotides in window %d\n", i_win);
			//getc(stdin);
		}
		else
		{
			// This is the average of poisson distribution for the heights in this window.
			double avg_read_depth = (double)(n_reads_per_window[i_win] * enrichment_mapped_fragment_length) / (double)cur_chr_data->n_uniquely_mappable_nucs_per_meg_win->at(i_win);
			double lambda = avg_read_depth;
			double log_lambda = xlog(lambda);

			double log_target_cdf = xlog(1.0 - target_fdr);

			// Update the log_cdf value for thresholds smaller than min_threshold.
			double current_log_cdf = -1.0 * lambda;
			double current_factor = -1.0 * lambda; // This is the probability value differentially updated in CDF computation.
			for(int thresh = 1; 
				thresh < min_thresh;
				thresh++)
			{
				double new_multiplier = xlog_div(log_lambda, xlog((double)thresh));
				current_factor = xlog_mul(current_factor, new_multiplier);
				current_log_cdf = xlog_sum(current_log_cdf, current_factor);
			}

			// For thresholds between min and max thresholds, compare the value of (log)CDF with the (log)target FDR.
			bool found_thresh = false;
			for(int thresh = min_thresh; 
				thresh <= max_thresh && !found_thresh; 
				thresh++)
			{
				double new_multiplier = xlog_div(log_lambda, xlog((double)thresh));
				current_factor = xlog_mul(current_factor, new_multiplier);
				current_log_cdf = xlog_sum(current_log_cdf, current_factor);

				if(current_log_cdf > log_target_cdf)
				{
					found_thresh = true;
					thresholds_per_win[i_win] = thresh; // Set the threshold in this window.
				}
			} // thresh loop.

			if(thresholds_per_win[i_win] > 0)
			{
				if(__DUMP_POISSON_BCKGRND_MSGS__)
					fprintf(stderr, "Window %d: %lf\n", i_win, thresholds_per_win[i_win]);
			}
		} // Check for positive number of reads in the window.
	} // i_win loop.

	delete [] n_reads_per_window;

	// Return the list of threshold arrays.
	return(thresholds_per_win);
}


