#include <stdio.h>
#include <stdlib.h>
#include "../../../lib/database/mapped_read/mapped_read_tools.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "simulated_background.h"
#include "poisson_background.h"
#include "../../../lib/chromosome/chromosome.h"
#include "compare_signal_tracks.h"
#include "../../../lib/database/mapped_read/fragment.h"
#include "peakseq.h"
#include "peak_region.h"
#include "chip_seq_chr_data.h"
#include "peakseq_output.h"
#include "../../../lib/utils/xmath/log/xlog_math.h"
#include <time.h>
#include <algorithm>
#include "../../../lib/utils/file/utils.h"

void peakseq(char* exp_id,
			char* chr_list_fp,
			int enrichment_fragment_length,
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

	int rand_seed = (simulation_seed != -1)?(simulation_seed):(time(0));
	srand(rand_seed);

	// Dump the seed file to regenerate the results.
	FILE* f_seed = open_f("seed.txt", "a");
	fprintf(f_seed, "%d\n", rand_seed);
	fclose(f_seed);
/*
    int enrichment_fragment_length = 200;
    double target_FDR = 0.05;
    int min_thresh = 1;
    int max_thresh = 100;
    int n_sims = 10;
    int min_gap_per_bw_peaks = 200;
    double P_f = 0.0; // Fraction of candidate peaks to exclude.
*/
    // Load the list of mappability map file paths.
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

/*
	clock_t beg = clock();
	vector<t_fragment*>* fore_frags = new vector<t_fragment*>();
	vector<t_fragment*>* rev_frags = new vector<t_fragment*>();

	load_fragments(false, "chr1.fa_mapped_reads.txt", fore_frags, rev_frags);
	//load_fragments_binary("chr1.fa_mapped_reads.bin", fore_frags, rev_frags);

	vector<t_fragment*>* chr_chip_seq_fragments = forwardize_combine_sort_fore_rev_strand_frags(fore_frags, rev_frags, 200);
	clock_t end = clock();
	printf("Loaded %d fragments. (%lf secs)\n", chr_chip_seq_fragments->size(), ((double)end - (double)beg) / CLOCKS_PER_SEC);

	//validate_dump_binary_files(chr_fps, false);

	exit(0);
*/

        // Read the ELAND file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
        //parse_ELAND_mapped_reads_file(chr_fps, chip_seq_eland_op_fp);

        // Read the ELAND file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
        //parse_ELAND_mapped_reads_file(chr_fps, control_eland_op_fp, true);

	//char mappability_map_fp[] = "Mapability_HG.txt";

	// Initialize the output object.
	//char op_fn[1000];

	//sprintf(op_fn, "%s_peaks", exp_id);


/*
    int peak_max_i = (int)floor(((double)cur_peak_region->max_height_start_i + (double)cur_peak_region->max_height_end_i) / 2);
    sprintf(op_msg, "%d %d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d",
            cur_peak_region->start_i_nuc,
            cur_peak_region->end_i_nuc,
            (int)cur_peak_region->n_chip_seq_reads,
            (int)cur_peak_region->n_adj_control_reads,
                n_excess,
                enrich,
            //xexp(cur_peak_region->log_p_val),
                cur_peak_region->log_p_val,
                cur_peak_region->q_val,
                peak_max_i,
                cur_peak_region->max_height);
*/

	// Allocate the final list of peaks, these are not necesssarily subset of union of all the peaks at all threshold.
	// Some of the peaks are cut and merged at window boundaries.
	vector<t_peak_region*>* annotated_peaks = new vector<t_peak_region*>();

	for(int i_chr = 0; i_chr < chr_fps->size(); i_chr++)
	{
		t_chip_seq_chr_data* cur_chr_data = new t_chip_seq_chr_data(chip_seq_reads_data_dirs, input_reads_data_dirs, chr_fps->at(i_chr), mappability_map_fp, enrichment_fragment_length);

		// Make sure that there are tags that are mapped on this chromosome.
		// Verify the input for this chromosome before starting the peak scoring.
		if(cur_chr_data->chip_seq_fragments->size() == 0 ||
			cur_chr_data->control_fragments->size() == 0 ||
			cur_chr_data->n_uniquely_mappable_nucs_per_meg_win->size() == 0)
		{
			printf("No (mappable) reads in the chromosome %s, skipping.\n", chr_fps->at(i_chr));
		}
		else
		{
			// Compute the scaling factor.
			double scaling_factor = get_input_scaling_factor(cur_chr_data->n_meg_wins(), cur_chr_data->chip_seq_fragments, cur_chr_data->control_fragments, P_f, enrichment_fragment_length);

			// Compare the signal tracks.
			if(scaling_factor > 0.0)
			{
				fprintf(stderr, "Computing all the peak regions over the chromosome for 1-100");
				vector<t_peak_region*>** peak_regions_per_threshold = cur_chr_data->chip_seq_enr_prof->get_peaks_per_thresholds(min_thresh, max_thresh, min_gap_per_bw_peaks);
				fprintf(stderr, "Computed all the peak regions over the chromosome for 1-100");

				// Count the peaks per threshold per each megbase window: This precounting makes the simulations much faster.
				int** n_peaks_per_thresh_per_window = count_peaks_per_thresh_per_window(cur_chr_data->n_meg_wins(),
																		peak_regions_per_threshold,
																		min_thresh,
																		max_thresh);

				printf("Processing %s\n", chr_fps->at(i_chr));

				// Simulation: For each chromosome, for each windows, for each threshold, for 40000 times, generate N_reads random reads.
				// Generate the enrichment profile for this window.
				double* thresholds_per_win = NULL;

				if(strcmp(background_model, "Simulated") == 0)
				{
					thresholds_per_win = simulate(cur_chr_data, 
													n_peaks_per_thresh_per_window,
													n_sims, 
													target_FDR, 
													enrichment_fragment_length, 
													min_gap_per_bw_peaks, 
													min_thresh, 
													max_thresh);
				}
				else if(strcmp(background_model, "Poisson") == 0)
				{
					thresholds_per_win = get_poisson_thresholds(cur_chr_data,
																target_FDR,
																enrichment_fragment_length,
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
																						enrichment_fragment_length,
																						min_gap_per_bw_peaks,
				 	              														thresholds_per_win,
																						peak_profs_dump_dir);

				// Add the new peaks to annotated peaks.
				for(int i_peak = 0; i_peak < cur_chr_annotated_peaks->size(); i_peak++)
				{
					annotated_peaks->push_back(cur_chr_annotated_peaks->at(i_peak));
				} // i_peak loop.

				//// Dump the peak profile for forward and reverse strand regions.
				//if(peak_profs_dump_dir != NULL && 
				//	strlen(peak_profs_dump_dir) > 0)
				//{
				//	vector<t_fragment*>* fore_chip_seq_fragments = forwardize_combine_sort_fore_rev_strand_frags(cur_chr_data->cur_chr_chip_seq_fore_frags, NULL, 30);
				//	t_enrichment_profile* fore_enr_prof = new t_enrichment_profile(fore_chip_seq_fragments, 30);
				//	fore_enr_prof->dump_peak_profiles_per_list(annotated_peaks, peak_profs_dump_dir, "peak_%d_fore_profile.txt");

				//	vector<t_fragment*>* rev_chip_seq_fragments = forwardize_combine_sort_fore_rev_strand_frags(NULL, cur_chr_data->cur_chr_chip_seq_rev_frags, 30);
				//	t_enrichment_profile* rev_enr_prof = new t_enrichment_profile(rev_chip_seq_fragments, 30);
				//	rev_enr_prof->dump_peak_profiles_per_list(annotated_peaks, peak_profs_dump_dir, "peak_%d_rev_profile.txt");
				//}

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
				printf("Skipping chromosome %s, scaling factor negative.\n", chr_fps->at(i_chr));
			}
		} // Check for # of mappable reads in the chromosome.

		// Free the data that is associated with the current chromosome.
		delete(cur_chr_data);
	} // i_chr loop.

	/*
	Following is the last step of the processing of the determined peaks: 
	1. Sort with respect to increasing p-values.
	2. Compute Q-values.
	3. Sort with respect to Q-values.
	4. Dump the peaks which have Q-value smaller than the input threshold for the Q-values.
	*/

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
	sprintf(op_title, "%-15s %-12s %-12s %-11s %-11s %-6s %-12s %-15s %-8s %-12s %-10s", "Chromosome", "Begin", "End", "#_Reads", "#_IP_Reads", "Excess", "Enrichment", "Log_p-value", "Q-value", "Max_Position", "Max_Height");
	t_peakseq_output::dump_PeakSeq_title(op_title);

	// Sort the peak regions with respect to their p-values, then compute Q-values: Note that after this point the peaks are not consecutive any more.
	sort(annotated_peaks->begin(), annotated_peaks->end(), sort_peak_region_per_log_p_val);

	printf("Processing %d annotated peaks.\n", annotated_peaks->size());

	for(int i_peak = 0; i_peak < annotated_peaks->size(); i_peak++)
	{
		double cur_peak_rank = (double)i_peak + 1;
		
		// Compute the Q-value for this annotated peak.
		annotated_peaks->at(i_peak)->q_val = (double)xexp(annotated_peaks->at(i_peak)->log_p_val) * ((double)annotated_peaks->size() / cur_peak_rank);
	} // i_peak 


	// Now sort the peaks with respect to Q-values.
	sort(annotated_peaks->begin(), annotated_peaks->end(), sort_peak_region_per_q_val);

	// Dump peak region information in extended PeakSeq output and in the narrowPeak output.
	for(int i_peak = 0; i_peak < annotated_peaks->size(); i_peak++)
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
			peak_max_i > cur_peak_region->end_i_nuc;
		}

/*
		if(cur_peak_region->max_height_start_i != cur_peak_region->max_height_end_i)
		{
			printf("peak<%d-%d> has a plateau\n", cur_peak_region->start_i_nuc, cur_peak_region->end_i_nuc);
			//getc(stdin);
		}
*/

/*
chr20           30150607     30151011     306         1           305    918.439237 -207.066085     0.000000 30150815     309       
chr20           31240682     31241252     305         1           304    1099.312480 -206.376190     0.000000 31240903     267
*/

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
					cur_peak_region->start_i_nuc,
					cur_peak_region->end_i_nuc,
					cur_peak_region->enrichment,
					-1.0 * cur_peak_region->log_p_val / xlog(10.0), // Convert log_2 into log_10.
					cur_peak_region->q_val,
					peak_max_i - cur_peak_region->start_i_nuc);

			t_peakseq_output::dump_narrowPeak_line(op_msg);
		}		
	} // i_peak.

	// Free memory for the peaks.
	for(int i_peak = 0; i_peak < annotated_peaks->size(); i_peak++)
	{
		delete(annotated_peaks->at(i_peak));
	} // i_peak loop.

	annotated_peaks->clear();
	delete(annotated_peaks);
}


