#include <stdio.h>
#include <stdlib.h>
#include "../../../lib/chromosome/chromosome.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "../../../lib/database/mapped_read/fragment.h"
#include <algorithm>
#include "simulated_background.h"
#include <math.h>
#include "compare_signal_tracks.h"
#include "peak_region.h"
#include "../../../lib/utils/xmath/log/xlog_math.h"
#include "chip_seq_chr_data.h"
#include "peakseq_output.h"
#include <string.h>


vector<t_peak_region*>* compare_signal_tracks(t_chip_seq_chr_data* cur_chr_data,
							vector<t_peak_region*>** peak_regions_per_threshold,
							int** n_peaks_per_thresh_per_window,
							double scaling_factor,
							int enrichment_fragment_length,
							int min_gap_per_bw_peaks,
							double* thresholds_per_win,
							char* peak_profs_dump_dir)
{
	//int i = 0;
	//double  n_control_reads = count_read_in_extended_region(56637247, 56637446, i, cur_chr_data->control_fragments, enrichment_fragment_length);
	//exit(0);

/*
	int i = 0;
        int n_control_reads = count_read_in_region(1235336, 1236193, i, cur_chr_control_fragments, enrichment_fragment_length);

	printf("%d-%d: %d control fragments in the extended region.\n", 1235336, 1236193, n_control_reads);
	exit(0);
*/

	int min_thresh = 1;
	int max_thresh = 100;
	printf("Comparing signal tracks!\n");

/*
	for(int i_peak = 0; i_peak < peak_regions_per_threshold[7 - min_thresh]->size(); i_peak++)
	{
		if(peak_regions_per_threshold[7 - min_thresh]->at(i_peak)->start_i_nuc == 36265272)
		{
			printf("Found the peak <%d-%d>\n", peak_regions_per_threshold[7 - min_thresh]->at(i_peak)->start_i_nuc, 
							peak_regions_per_threshold[7 - min_thresh]->at(i_peak)->end_i_nuc);
			
		}
	}

	cur_chr_data->chip_seq_enr_prof->dump_profile("chrX_profile.txt", 36 * MEG_BASE, 36 * MEG_BASE + 500000);
*/

	// Initialize all the indices to 0.
	int* last_peak_list_indices = new int[max_thresh - min_thresh + 1];
	for(int i_thr = min_thresh; i_thr <= max_thresh; i_thr++)
	{
		last_peak_list_indices[i_thr - min_thresh] = 0;
	}

	vector<t_peak_region*>* cur_chr_annotated_peaks = new vector<t_peak_region*>();

	// For each window, use the thresholds and determine the peaks, using the number of reads in the identified peaks, based on the number 
	// of reads, determine significance of the peak.
	int n_wins = cur_chr_data->n_meg_wins();
	int last_peak_end = -10000;
	int last_peak_start = 0;

	int last_chip_seq_frag_i = 0;
	int last_control_frag_i = 0;
	for(int i_win = 0; i_win < n_wins; i_win++)
	{
		// Get the peaks in this window: Use the threhold for this window. It is necessary to account for the first peak. 
		int cur_win_thresh = (int)thresholds_per_win[i_win];

		if(cur_win_thresh > 0 &&
			n_peaks_per_thresh_per_window[cur_win_thresh - min_thresh][i_win] > 0)
		{
			process_peaks_per_window(i_win,
					cur_win_thresh,
					cur_chr_data,
					peak_regions_per_threshold[cur_win_thresh - min_thresh],
                    last_peak_list_indices[cur_win_thresh - min_thresh],
                    last_peak_start,
                    last_peak_end,
                    min_gap_per_bw_peaks,
					scaling_factor,
					enrichment_fragment_length,
					last_chip_seq_frag_i,
					last_control_frag_i,
					cur_chr_annotated_peaks);
		}
		else
		{
		}

/*
		if(i_win == 56)
		{
			getc(stdin);
		}
*/
		// Determine the list of peaks in this window.
	} // i_win loop.

	int i_last_win = 0;

    for(int i_win = 0; i_win < n_wins; i_win++)
    {
        // Get the peaks in this window: Use the threhold for this window. It is necessary to account for the first peak.
        int cur_win_thresh = (int)thresholds_per_win[i_win];

		if(cur_win_thresh > 0)
		{
			if(n_peaks_per_thresh_per_window[cur_win_thresh - min_thresh][i_win] > 0)
			{
				// For this window, if there are peaks, then this is the last window.
				i_last_win = i_win;
			}
		}
        else
        {
        }
        // Determine the list of peaks in this window.
    } // i_win loop.

	printf("Last peak is <%d, %d>\n", peak_regions_per_threshold[(int)thresholds_per_win[i_last_win] - min_thresh]->back()->start_i_nuc, peak_regions_per_threshold[(int)thresholds_per_win[i_last_win] - min_thresh]->back()->end_i_nuc);

	// Dump the last peak of the last window that has peaks in it.
    // This ensures that the very last peak is dumped.
	t_peak_region* new_annot_peak = dump_region_info(cur_chr_data->chr_id,
													i_last_win,
													//peak_regions_per_threshold[(int)thresholds_per_win[i_last_win] - min_thresh]->back()->start_i_nuc,
													//peak_regions_per_threshold[(int)thresholds_per_win[i_last_win] - min_thresh]->back()->end_i_nuc,
													last_peak_start,
													last_peak_end,
													cur_chr_data->chip_seq_fragments,
													cur_chr_data->control_fragments,
													scaling_factor,
													enrichment_fragment_length,
													last_chip_seq_frag_i,
    												last_control_frag_i);

	cur_chr_annotated_peaks->push_back(new_annot_peak);

    // Update the maximums for the annotated peaks.
    cur_chr_data->chip_seq_enr_prof->set_max_per_peak_region_list(cur_chr_annotated_peaks);

 //       // Dump the peak profile.
	//if(peak_profs_dump_dir != NULL && 
	//	strlen(peak_profs_dump_dir) > 0)
	//{
	//        cur_chr_data->chip_seq_enr_prof->dump_peak_profiles_per_list(annotated_peaks, peak_profs_dump_dir, "peak_%d_%d_chip_seq_profile.txt");
	//}

	return(cur_chr_annotated_peaks);
}

bool sort_peak_region_per_log_p_val(t_peak_region* peak1, t_peak_region* peak2)
{
	// Smaller p-value means more significant.
	return(peak1->log_p_val < peak2->log_p_val);
}

bool sort_peak_region_per_q_val(t_peak_region* peak1, t_peak_region* peak2)
{
	if(peak1->q_val < peak2->q_val)
	{
		return(true);
	}
	else if(peak1->q_val == peak2->q_val)
	{
		return(peak1->log_p_val < peak2->log_p_val);
	}
	else
	{
		return(false);
	}
}

void process_peaks_per_window(int i_win, 
				int cur_win_threshold, 
				t_chip_seq_chr_data* cur_chr_data,
				vector<t_peak_region*>* thr_peaks,
				int& last_peak_list_index, // This is the peak index for the current threshold to decrease the search time of the peaks in the window.
				int& last_peak_start, // The starting index for the last identified peak
				int& last_peak_end, // The end index for the last identified peak.
				int min_gap_per_bw_peaks,
				double scaling_factor,
				int enrichment_fragment_length,
				int& last_chip_seq_frag_i, // Global last index for the ChIP-Seq enriched fragments, this decreases the search time when the fragments are counted in the peak identified in the chip-seq signal. Not updated in the extended searches.
				int& last_control_frag_i,
				vector<t_peak_region*>* annotated_peaks) // Global last index for the control enriched fragments for efficient searching.
{	
	// Find ihe next peak that is overlapping with this window.
	fprintf(stderr, "Processing the %d peaks for window %d with threshold %d, last peak was <%d-%d>\n", thr_peaks->size(), i_win, cur_win_threshold, last_peak_start, last_peak_end);
	//getc(stdin);

	printf("last_peak_list_index is %d\n", last_peak_list_index);

	//printf("End nucleotide of first peak is %d\n", thr_peaks->at(last_peak_list_index)->end_i_nuc);

	// Find the first peak in the window, look for the overlap between beginning part of the window and end of the peak.
	while(last_peak_list_index < thr_peaks->size() && 
		(thr_peaks->at(last_peak_list_index)->end_i_nuc < (i_win * MEG_BASE)))
	{
		//printf("Skipping %d. peak ending at %d.\n", last_peak_list_index, thr_peaks->at(last_peak_list_index)->end_i_nuc);
		last_peak_list_index++;
	} // peak search loop.

	if(last_peak_list_index == thr_peaks->size())
	{
		return;
	}

	// At this point, there are 5 possibilities for the peak's positioning with respect to the window and 4 of them must be handled. 
	// It may be a good idea to inspect if the following code handles all of them properly.

	// If the peaks list is finihsed or this peak does not overlap with this window, do not process the peaks in this window.
	if(last_peak_list_index < thr_peaks->size() &&
		(thr_peaks->at(last_peak_list_index)->start_i_nuc >= (i_win+1) * MEG_BASE))
	{
		printf("There are no peaks in the window %d for the threshold %d.\n", i_win, cur_win_threshold);

		if(last_peak_list_index > 0)
		{
			last_peak_list_index--;
		}
		return;
		//exit(0);
	}

	// Determine the indices for the overlapping parts of this peak with the current window.
	int first_peak_win_start = (thr_peaks->at(last_peak_list_index)->start_i_nuc < (i_win * MEG_BASE)) ? (i_win * MEG_BASE) : (thr_peaks->at(last_peak_list_index)->start_i_nuc );
	int first_peak_win_end = (thr_peaks->at(last_peak_list_index)->end_i_nuc >= ((i_win+1) * MEG_BASE)) ? ((i_win+1) * MEG_BASE-1) : (thr_peaks->at(last_peak_list_index)->end_i_nuc);

	// Is this peak far enough to be a peak by itself?
	if(last_peak_list_index < thr_peaks->size() &&
		first_peak_win_start > (last_peak_end + min_gap_per_bw_peaks))
	{
		if(last_peak_end >= 0)
		{
			// Dump the last peak.
			printf("Dumping the information for peak<%d-%d>\n", last_peak_start, last_peak_end);
	
			t_peak_region* new_annot_peak = dump_region_info(cur_chr_data->chr_id,
															i_win,
															last_peak_start,
															last_peak_end,
															cur_chr_data->chip_seq_fragments,
															cur_chr_data->control_fragments,
															scaling_factor,
															enrichment_fragment_length,
															last_chip_seq_frag_i,
															last_control_frag_i);

			// Add this peak to annotated peak list.
			annotated_peaks->push_back(new_annot_peak);
		}

		last_peak_start = first_peak_win_start; // Reset the start for the last peak.
	}
	else
	{
		// Must merge the previous peak in the previous window.
	}

	// Set the end of the peak to the end of the currently found peak, the first peak in the window.
	// There are two possibilities:
	// 1. The last peak in the previous window is far enough and is dumped already in the above conditional.
	// 2. The last peak was not far enough and the above conditional was not entered. In this case the following line "merges" the peaks by extending the end of last read peak, thereby 
	// when the next peak is dumped, it is dumped with the extended indices. This may introduce artifical discontinuities. It may be a good idea to ignore these boundary peaks while processing
	// the peaks.
	last_peak_end = first_peak_win_end;
	last_peak_list_index++;	// Move to the next peak.

	// Now go over all the remaining peaks that have overlapping parts with this window: Note that the min gap check is not required for the identified peaks.
	while(last_peak_list_index < thr_peaks->size() &&
		(thr_peaks->at(last_peak_list_index)->start_i_nuc < (i_win+1) * MEG_BASE))
	{
		// Dump the information for the last identified region. This region may be formed by merging parts of two peaks.
                //printf("Dumping the information for peak<%d-%d>\n", last_peak_start, last_peak_end);
		t_peak_region* new_annot_peak = dump_region_info(cur_chr_data->chr_id, i_win,
														last_peak_start,
														last_peak_end,	
														cur_chr_data->chip_seq_fragments,
														cur_chr_data->control_fragments,
														scaling_factor,
														enrichment_fragment_length,
														last_chip_seq_frag_i,
														last_control_frag_i);
		// Add this peak to annotated peak list.
		annotated_peaks->push_back(new_annot_peak);

		last_peak_start = thr_peaks->at(last_peak_list_index)->start_i_nuc;
		last_peak_end = (thr_peaks->at(last_peak_list_index)->end_i_nuc >= ((i_win+1) * MEG_BASE)) ? ((i_win+1) * MEG_BASE-1) : (thr_peaks->at(last_peak_list_index)->end_i_nuc);

		// Move to the next peak for this threshold.
		last_peak_list_index++;
	} // Go over all the peaks in the window until the peaks are finished or the peaks do not overlap with the window.

	// Make sure the peak index for this threshold 
	if(last_peak_list_index > 0)
	{
		last_peak_list_index--;
	}

	// Note that at this point, the last peak is not dumped, yet. The reason for this is simple: The peak may need to be merged with the first peak in the next window.
	// However, the last peak of the last window has to be manually dumped at the end of the window loop.

	//exit(0);	
}

t_peak_region* dump_region_info(char* chr_id,
								int i_win, 
								int start, 
								int end,
								vector<t_fragment*>* chip_seq_fragments, 
								vector<t_fragment*>* control_fragments, 
								double scaling_factor,
								int enrichment_fragment_length,
								int& last_chip_seq_frag_i,
								int& last_control_frag_i)
{
	printf("Dumping information for region<%d-%d>\n", start, end);

	// Number of reads in the region.
	int n_chip_seq_reads = 0;
	double n_control_reads = 0;

	int rep_cnt = 0;

	n_chip_seq_reads = count_read_in_region(start, end, last_chip_seq_frag_i, chip_seq_fragments, enrichment_fragment_length);
	
	if(end - start > MAX_REG_EXT)
	{
		n_control_reads = (double)count_read_in_region(start, end, last_control_frag_i, control_fragments, enrichment_fragment_length);
	}
	else
	{
		n_control_reads = count_read_in_extended_region(start, end, last_control_frag_i, control_fragments, enrichment_fragment_length);
	}

	double adjusted_n_control_reads = ceil(n_control_reads * scaling_factor);

        if (n_control_reads == 0.0)
                n_control_reads = 1.0;
        if (adjusted_n_control_reads == 0.0)
                adjusted_n_control_reads = 1.0;

	int n = (int)adjusted_n_control_reads + n_chip_seq_reads;

	double* logs = new double[n + 5];
	buffer_logs(logs, n+2);

	// Logarithm of 1/2.
	double log_half = xlog(0.5);

	// Compute the null hypothesis probability.
	double log_p_val = xlog(0.0);

	// Compute the factroial of n since it will be needed constantly.
	double* log_fact_buffer = new double[n + 5];
	buffer_log_facts(log_fact_buffer, n+2, logs);

	// Compute each possibility.
	for(int i = 0; i <= (int)adjusted_n_control_reads; i++)
	{
		// Use trick for computing power.
		double log_cur_half_pow = log_half * n;
		double log_cur_perm = 0.0; // = xlog(1.0).

		// Compute the current permutation, using shortcuts, bypassing the function calls from xlog_math library.
		log_cur_perm = log_fact_buffer[n] - (log_fact_buffer[i] + log_fact_buffer[n-i]);				

		log_p_val = xlog_sum(log_p_val, (log_cur_perm + log_cur_half_pow));
	} // i loop.

	// Free arrays.
	delete[] logs;
	delete[] log_fact_buffer;
    printf("%d sample tags, %lf control tags, %lf adjusted control tags, pval = %lf (%lf)\n", n_chip_seq_reads, n_control_reads, adjusted_n_control_reads, xexp(log_p_val), log_p_val);

	// Allocate the peak region and set it.
	t_peak_region* new_peak_region = new t_peak_region(start, end);
	new_peak_region->log_p_val = log_p_val;
    new_peak_region->q_val = 0.0;
    new_peak_region->n_chip_seq_reads = n_chip_seq_reads;
	new_peak_region->n_control_reads = n_control_reads;
	new_peak_region->n_adj_control_reads = adjusted_n_control_reads;
	new_peak_region->enrichment = (double)n_chip_seq_reads / ((double)n_control_reads * scaling_factor);
	new_peak_region->chr_id = new char[strlen(chr_id) + 2];
	strcpy(new_peak_region->chr_id, chr_id);
	//new_peak_region->excess = new_peak_region->n_chip_seq_reads - new_peak_region->n_adj_control_reads;
	//new_peak_region->

/*
    char op_msg[500];
    //sprintf(op_msg, "%d %d\t%d\t%d\t%lf (%lf)", start, end, n_chip_seq_reads, (int)adjusted_n_control_reads, xexp(log_p_val), log_p_val);
    sprintf(op_msg, "%d %d\t%d\t%d\t%lf",
            new_peak_region->start_i_nuc,
            new_peak_region->end_i_nuc,
            (int)new_peak_region->n_chip_seq_reads,
            (int)new_peak_region->n_adj_control_reads,
            xexp(new_peak_region->log_p_val));

    t_peakseq_output::dump_PeakSeq_line(op_msg);
*/	
	return(new_peak_region);
}

int count_read_in_region(int start, int end, int& last_frag_list_i, vector<t_fragment*>* frag_list, int enrichment_fragment_length)
{
	int n_region_reads = 0;
        int last_added_index = -1000;
        int rep_cnt = 0;

        for(int i_frag = last_frag_list_i; i_frag < frag_list->size(); i_frag++)
        {
                if(frag_list->at(i_frag)->base_index > end)
                {
                        break;
                }
                else if(frag_list->at(i_frag)->base_index + enrichment_fragment_length >= start)
                {
			// Update the last frag list index.
			last_frag_list_i = i_frag;

                        // Update repeat counter for this fragment, if necessary.
                        if(last_added_index == frag_list->at(i_frag)->base_index)
                        {
                                rep_cnt++;
                        }
                        else
                        {
                                last_added_index = frag_list->at(i_frag)->base_index;
                                rep_cnt = 1;
                        }

                        // If repeat counter is greater than 3, do not include this fragment.
                        if(rep_cnt > 3)
                        {
                        }
                        else
                        {
                                n_region_reads++;
                        }
                }
        } // i_ch_frag loop.

	return(n_region_reads);
}

double count_read_in_extended_region(int start, int end, int& last_frag_list_i, vector<t_fragment*>* frag_list, int enrichment_fragment_length)
{
	printf("Counting extended region tags: %d\n", last_frag_list_i);

	int cur_last_frag_list_i = last_frag_list_i;
	int n_region_reads_int = count_read_in_region(start, end, last_frag_list_i, frag_list, enrichment_fragment_length);

	int ext_start = (start > EXTENDED_REGION_SIZE) ? start - EXTENDED_REGION_SIZE : 1;
	int ext_end = end + EXTENDED_REGION_SIZE;

	printf("Before extended region count index is %d\n", cur_last_frag_list_i);

	// Go back and find the first index that is before the ext_start. Note that there may be multiple fragments that have same position, 
	// and we want to count at most 3 of those, so the check here does GEQ, not just GT.
	while((frag_list->at(cur_last_frag_list_i)->base_index + enrichment_fragment_length) >= ext_start &&
		cur_last_frag_list_i > 0)
	{
		cur_last_frag_list_i--;
	}

	int n_ext_region_reads_int = count_read_in_region(ext_start, ext_end, cur_last_frag_list_i, frag_list, enrichment_fragment_length);

	double n_ext_region_reads_ld = (double)n_ext_region_reads_int * ( ((double) (end - start + 1)) /
							(end - start + 1 + (2 * EXTENDED_REGION_SIZE)) );

        printf("Extended region count: %d x %lf = %lf\n", n_ext_region_reads_int,  ( ((double) (end - start + 1)) /
                                                        (end - start + 1 + (2 * EXTENDED_REGION_SIZE) )), n_ext_region_reads_ld);

	// Return the maximum of the adjusted ext_control_tags and the normal
	// control_tags.
	double n_region_reads_ld = (n_ext_region_reads_ld > (double)n_region_reads_int) ? n_ext_region_reads_ld : (double)n_region_reads_int;
	return(n_region_reads_ld);
}

void buffer_logs(double* buff, int n)
{
	for(int i = 0; i <= n; i++)
	{
		buff[i] = xlog(i);
	} // i loop.
}

void buffer_log_facts(double* log_fact_buffer, int n, double* logs_buffer)
{
	log_fact_buffer[0] = 0.0f; // 0! = 1.
	for(int i = 1; i <= n; i++)
	{
		// Bypass xlog_math functions.
		log_fact_buffer[i] = logs_buffer[i] + log_fact_buffer[i-1];
	}
}


