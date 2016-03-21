#ifndef __COMPARE_SIGNAL_TRACKS__
#define __COMPARE_SIGNAL_TRASKS__

#include <vector>
using namespace std;

#define MAX_REG_EXT (1000)
#define EXTENDED_REGION_SIZE (2000)

class t_chromosome;
struct t_mapped_fragment;
class t_enrichment_profile;

int count_read_in_region(int start, int end, int& last_frag_i, vector<t_mapped_fragment*>* frag_list, int enrichment_mapped_fragment_length);
double count_read_in_extended_region(int start, int end, int& last_frag_i, vector<t_mapped_fragment*>* frag_list, int enrichment_mapped_fragment_length);

void buffer_logs(double* buff, int n);
void buffer_log_facts(double* log_fact_buffer, int n, double* logs_buffer);

bool sort_peak_region_per_log_p_val(t_peak_region* peak1, t_peak_region* peak2);
bool sort_peak_region_per_q_val(t_peak_region* peak1, t_peak_region* peak2);

vector<t_peak_region*>* compare_signal_tracks(t_chip_seq_chr_data* cur_chr_data,
													//vector<t_peak_region*>* annotated_peaks,
													vector<t_peak_region*>** peak_regions_per_threshold,
													int** n_peaks_per_thresh_per_window,
													double scaling_factor,
													int enrichment_mapped_fragment_length,
													int min_gap_per_bw_peaks,
													double* thresholds_per_win,
													char* peak_profs_dump_dir,
													int min_thresh,
													int max_thresh);

void process_peaks_per_window(int i_win,
                                int cur_win_threshold,
				t_chip_seq_chr_data* cur_chr_data,
                                vector<t_peak_region*>* thr_peaks,
                                int& last_peak_list_index,
                                int& last_peak_start,
                                int& last_peak_end,
                                int min_gap_per_bw_peaks,
                                double scaling_factor,
				int enrichment_mapped_fragment_length,
                                int& last_chip_seq_frag_i,
                                int& last_control_frag_i,
				vector<t_peak_region*>* annotated_peaks);


t_peak_region* dump_region_info(char* chr_id, 
								int i_win,
	                        int start,
	                        int end,
	                        vector<t_mapped_fragment*>* chip_seq_fragments,
	                        vector<t_mapped_fragment*>* control_fragments,
	                        double scaling_factor, 
				int enrichment_mapped_fragment_length,
				int& last_chip_seq_frag_i,
				int& last_control_frag_i);

#endif // __COMPARE_SIGNAL_TRASKS__



