#include <stdio.h>
#include <stdlib.h>
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "simulated_background.h"
#include "poisson_background.h"
#include "genomics_coords.h"
#include "compare_signal_tracks.h"
#include "peakseq.h"
#include "peak_region.h"
#include "chip_seq_chr_data.h"
#include "peakseq_output.h"
#include "xlog_math.h"
#include <time.h>
#include <algorithm>
#include <string.h>
#include "utils.h"

bool __DUMP_PEAKSEQ_MSGS__ = false;

void peakseq(char* exp_id,
			char* chr_list_fp,
			int enrichment_mapped_fragment_length,
			char* background_model,
			double target_FDR, // FDR threshold for initial height threshold selection.
			double max_Qval, // Q-value threshold.
			int min_thresh,
			int max_thresh,
			int n_sims,
			int min_gap_per_bw_peaks,
			double P_f,
			vector<char*>* chip_seq_reads_data_dirs,
			vector<char*>* input_reads_data_dirs,
			char* mappability_map_fp,
			char* peak_profs_dump_dir,
			char* narrowPeak_op_fp,
			int simulation_seed)
{
	int rand_seed = (simulation_seed != -1)?(simulation_seed):((int)time(0));
	srand(rand_seed);

	// Dump the seed file to regenerate the results.
	FILE* f_seed = open_f("seed.txt", "a");
	fprintf(f_seed, "%d\n", rand_seed);
	fclose(f_seed);

    // Load the chromosome id's to process.
    vector<char*>* chr_fps = new vector<char*>();
    FILE* f_chr_list = open_f(chr_list_fp, "r");
    if(f_chr_list == NULL)
    {
            printf("Could not open chrososome list file %s\n", chr_list_fp);
            exit(0);
    }

    char cur_fp[100];
    while(fscanf(f_chr_list, "%s", cur_fp) == 1)
    {
            char* new_chr_fp = new char[strlen(cur_fp) + 2];
            strcpy(new_chr_fp, cur_fp);
            chr_fps->push_back(new_chr_fp);
            //printf("%s\n", new_chr_fp);
    }
    fclose(f_chr_list);

	// Allocate the final list of peaks, these are not necesssarily subset of union of all the peaks at all threshold.
	// Some of the peaks are cut and merged at window boundaries.
	vector<t_peak_region*>* annotated_peaks = new vector<t_peak_region*>();

	for(int i_chr = 0; i_chr < (int)chr_fps->size(); i_chr++)
	{
		fprintf(stderr, "%s..", chr_fps->at(i_chr));

		t_chip_seq_chr_data* cur_chr_data = new t_chip_seq_chr_data(chip_seq_reads_data_dirs, input_reads_data_dirs, chr_fps->at(i_chr), mappability_map_fp, enrichment_mapped_fragment_length);

		// Make sure that there are tags that are mapped on this chromosome.
		// Verify the input for this chromosome before starting the peak scoring.
		if(cur_chr_data->chip_seq_fragments->size() == 0 ||
			cur_chr_data->control_fragments->size() == 0 ||
			cur_chr_data->n_uniquely_mappable_nucs_per_meg_win->size() == 0)
		{
if(__DUMP_PEAKSEQ_MSGS__)
			fprintf(stderr, "No (mappable) reads in the chromosome %s, skipping.\n", chr_fps->at(i_chr));
		}
		else
		{
			// Compute the scaling factor.
			double scaling_factor = get_input_scaling_factor(cur_chr_data->n_meg_wins(), cur_chr_data->chip_seq_fragments, cur_chr_data->control_fragments, P_f, enrichment_mapped_fragment_length);

			// Compare the signal tracks.
			if(scaling_factor > 0.0)
			{
if(__DUMP_PEAKSEQ_MSGS__)
				fprintf(stderr, "Computing all the peak regions over the chromosome for 1-100");

				vector<t_peak_region*>** peak_regions_per_threshold = cur_chr_data->chip_seq_enr_prof->get_peaks_per_thresholds(min_thresh, max_thresh, min_gap_per_bw_peaks);

if(__DUMP_PEAKSEQ_MSGS__)
				fprintf(stderr, "Computed all the peak regions over the chromosome for 1-100");

				// Count the peaks per threshold per each megbase window: This precounting makes the simulations much faster.
				int** n_peaks_per_thresh_per_window = count_peaks_per_thresh_per_window(cur_chr_data->n_meg_wins(),
																		peak_regions_per_threshold,
																		min_thresh,
																		max_thresh);			

				// Simulation: For each chromosome, for each windows, for each threshold, for 40000 times, generate N_reads random reads.
				// Generate the enrichment profile for this window.
				double* thresholds_per_win = NULL;

				if(strcmp(background_model, "Simulated") == 0)
				{
					thresholds_per_win = simulate(cur_chr_data, 
													n_peaks_per_thresh_per_window,
													n_sims, 
													target_FDR, 
													enrichment_mapped_fragment_length, 
													min_gap_per_bw_peaks, 
													min_thresh, 
													max_thresh);
				}
				else if(strcmp(background_model, "Poisson") == 0)
				{
					thresholds_per_win = get_poisson_thresholds(cur_chr_data,
																target_FDR,
																enrichment_mapped_fragment_length,
																min_thresh,
																max_thresh);
				}
				else
				{
					printf("Could not resolve background model string %s\n", background_model);
					exit(0);
				}

				// Copare the signal tracks and compute significance for each peak.
				vector<t_peak_region*>* cur_chr_annotated_peaks = compare_signal_tracks(cur_chr_data,
																						peak_regions_per_threshold, 
																						n_peaks_per_thresh_per_window,
																						scaling_factor,
																						enrichment_mapped_fragment_length,
																						min_gap_per_bw_peaks,
				 	              														thresholds_per_win,
																						peak_profs_dump_dir,
																						min_thresh,
																						max_thresh);

				// Add the new peaks to annotated peaks.
				for(int i_peak = 0; i_peak < (int)cur_chr_annotated_peaks->size(); i_peak++)
				{
					annotated_peaks->push_back(cur_chr_annotated_peaks->at(i_peak));
				} // i_peak loop.

				// Free memory, note that it is possible that this chromosome is skipped.
				delete [] thresholds_per_win;
				free_n_peaks_per_thresh_per_window(n_peaks_per_thresh_per_window,
												min_thresh,
												max_thresh);

				free_peaks_per_threshold(peak_regions_per_threshold,
											min_thresh,
											max_thresh);
			} // positive scaling_factor check
			else
			{
				fprintf(stderr, "Skipping chromosome %s, scaling factor negative.\n", chr_fps->at(i_chr));
			}
		} // Check for # of mappable reads in the chromosome.

		// Free the data that is associated with the current chromosome.
		delete(cur_chr_data);
	} // i_chr loop.

	fprintf(stderr, "Done.\n");

	/*/
	Following is the last step of the processing of the determined peaks: 
	1. Sort with respect to increasing p-values.
	2. Compute Q-values.
	3. Sort with respect to Q-values.
	4. Dump the peaks which have Q-value smaller than the input threshold for the Q-values.
	/*/
	t_peakseq_output* ps_op = NULL;
	if(strlen(narrowPeak_op_fp) > 0)
	{
		ps_op = new t_peakseq_output(exp_id, narrowPeak_op_fp);
	}
	else
	{
		ps_op = new t_peakseq_output(exp_id);
	}

	char op_title[1000];

	// sprintf(op_msg, "%-12d %-12d %-11d %-11d %-6d %-12.3lf %-15.3lf %-8.3lf %-12d %-10d",
	sprintf(op_title, "%-15s %-12s %-12s %-11s %-11s %-6s %-12s %-15s %-8s %-12s %-10s", "Chromosome", "Begin", "End", "#_ChIP_Reads", "#_Control_Reads", "Excess", "Enrichment", "Log_p-value", "Q-value", "Max_Position", "Max_Height");
	t_peakseq_output::dump_PeakSeq_title(op_title);

	// Sort the peak regions with respect to their p-values, then compute Q-values: Note that after this point the peaks are not consecutive any more.
	sort(annotated_peaks->begin(), annotated_peaks->end(), sort_peak_region_per_log_p_val);

if(__DUMP_PEAKSEQ_MSGS__)
	printf("Processing %d annotated peaks.\n", (int)annotated_peaks->size());

	for(int i_peak = 0; i_peak < (int)annotated_peaks->size(); i_peak++)
	{
		double cur_peak_rank = (double)i_peak + 1;
		
		// Compute the Q-value for this annotated peak.
		annotated_peaks->at(i_peak)->q_val = (double)xexp(annotated_peaks->at(i_peak)->log_p_val) * ((double)annotated_peaks->size() / cur_peak_rank);
	} // i_peak 


	// Now sort the peaks with respect to Q-values.
	sort(annotated_peaks->begin(), annotated_peaks->end(), sort_peak_region_per_q_val);

	// Dump peak region information in extended PeakSeq output and in the narrowPeak output.
	for(int i_peak = 0; i_peak < (int)annotated_peaks->size(); i_peak++)
	{
		t_peak_region* cur_peak_region = annotated_peaks->at(i_peak);
	    char op_msg[1500];
	    //sprintf(op_msg, "%d %d\t%d\t%d\t%lf (%lf)", start, end, n_chip_seq_reads, (int)adjusted_n_control_reads, xexp(log_p_val), log_p_val);

		int n_excess = (int)(cur_peak_region->n_chip_seq_reads - cur_peak_region->n_adj_control_reads);
		//double enrich = (double)cur_peak_region->n_chip_seq_reads / ((double)cur_peak_region->n_control_reads * scaling_factor);

		int peak_max_i = (int)floor(((double)cur_peak_region->max_height_start_i + (double)cur_peak_region->max_height_end_i) / 2);

		if(peak_max_i > cur_peak_region->end_i_nuc ||
			peak_max_i > cur_peak_region->start_i_nuc)
		{
			peak_max_i = cur_peak_region->end_i_nuc;
		}

		// Check the Q-value threshold before dumping.
		if(cur_peak_region->q_val < max_Qval)
		{
			sprintf(op_msg, "%-15s %-12d %-12d %-11d %-11d %-6d %-12.3lf %-15.3lf %-8.3lf %-12d %-10d",
					cur_peak_region->chr_id,
					cur_peak_region->start_i_nuc,
					cur_peak_region->end_i_nuc,
					(int)cur_peak_region->n_chip_seq_reads,	
					(int)cur_peak_region->n_adj_control_reads,
				n_excess,
				cur_peak_region->enrichment,
					//xexp(cur_peak_region->log_p_val),
				cur_peak_region->log_p_val,
				cur_peak_region->q_val,
				peak_max_i - cur_peak_region->start_i_nuc, // Dump the relative distance of the peak max.
				cur_peak_region->max_height);

			t_peakseq_output::dump_PeakSeq_line(op_msg);

			// Dump the narrowPeak output line.
			sprintf(op_msg, "%-15s\t%-12d\t%-12d\t.\t0\t.\t%-12.3lf\t%e\t%e\t%d",
					cur_peak_region->chr_id,
					cur_peak_region->start_i_nuc - CODEBASE_COORDS::start_base + BED_COORDS::start_base,
					cur_peak_region->end_i_nuc - CODEBASE_COORDS::end_base + BED_COORDS::end_base,
					cur_peak_region->enrichment,
					-1.0 * cur_peak_region->log_p_val / xlog(10.0), // Convert log_2 into log_10.
					cur_peak_region->q_val,
					peak_max_i - cur_peak_region->start_i_nuc);

			t_peakseq_output::dump_narrowPeak_line(op_msg);
		}		
	} // i_peak.

	// Free memory for the peaks.
	for(int i_peak = 0; i_peak < (int)annotated_peaks->size(); i_peak++)
	{
		delete(annotated_peaks->at(i_peak));
	} // i_peak loop.

	annotated_peaks->clear();
	delete(annotated_peaks);
}


