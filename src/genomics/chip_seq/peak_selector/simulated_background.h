#ifndef __SIMULATE__
#define __SIMULATE__

#include <vector>
using namespace std;

class t_chromosome;
class t_enrichment_profile;
struct t_fragment;
class t_chip_seq_chr_data;

double* simulate(t_chip_seq_chr_data* cur_chr_data,
		int** n_peaks_per_thresh_per_window,
                int n_iterations,
                double target_fdr,
                int enrichment_fragment_length,
                int min_gap_per_bw_peaks,
                int min_thresh,
                int max_thresh);

#endif // __SIMULATE__
