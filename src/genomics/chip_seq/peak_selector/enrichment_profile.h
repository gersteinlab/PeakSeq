#ifndef __ENRICHMENT_PROFILE__
#define __ENRICHMENT_PROFILE__

#include <vector>
using namespace std;

class t_chromosome;
struct t_fragment;
class t_peak_region;
class t_chip_seq_chr_data;
struct t_annot_region;
//struct t_profile_site;
struct t_frag_end;

class t_enrichment_profile
{
public:
    // Segment length is the global length of the sequence over which the enrichment is performed, which may be experimental or computational.
    t_enrichment_profile(char* chr_id, vector<t_fragment*>* fragments, int enrichment_fragment_length);
    t_enrichment_profile(char* chr_id, t_fragment** fragments, int enrichment_fragment_length);
	t_enrichment_profile(char* chr_id, char* sgr_fp, int enrichment_fragment_length);
	t_enrichment_profile();

    ~t_enrichment_profile();

	char* chr_id;

	// A profile is stored as the fragment end points. Do not really need to store where the fragments are.
	vector<t_frag_end*>* frag_ends;

	void add_fragment_ends(int start, int length);

	// Count the peaks for this enrichment factor with the threshold.
	//int count_peaks(int threshold, int min_gap_bw_peaks, int min_i_nuc, int max_i_nuc);
	void count_peaks(int* n_peaks, int min_thresh, int max_thresh, int min_gap_bw_peaks, int min_i_nuc, int max_i_nuc);
	void count_peaks_per_frags(int* n_peaks_per_thresh, vector<int>* frag_begins, vector<int>* frag_ends, int min_thresh, int max_thresh, int min_gap_per_bw_peaks);

	vector<t_peak_region*>** get_peaks_per_thresholds(int min_threshold, int max_threshold, int min_gap_bw_peaks);

	void dump_profile(char* fp,
						int min_i_nuc,
						int max_i_nuc);

	//vector<t_profile_site*>* buffer_profile(int min_i_nuc, int max_i_nuc);

	// Dump the profiles for the list of peaks. This is more efficient than dumping each peak separately.
	void dump_peak_profiles_per_list(vector<t_peak_region*>* peak_list, char* dump_dir, char* file_name_fmt_str);

	// Set the maximum and the position of maximum for all the peak regions in the list using the enrichment profile.
	void set_max_per_peak_region_list(vector<t_peak_region*>* peak_region_list);
};

int** count_peaks_per_thresh_per_window(int n_meg_wins,
                                                vector<t_peak_region*>** peaks_per_threshold,
                                                int min_thresh,
                                                int max_thresh);

void free_n_peaks_per_thresh_per_window(int** n_peaks_per_thresh_per_window,
                                        int min_thresh,
                                        int max_thresh);


void free_peaks_per_threshold(vector<t_peak_region*>** peak_regions_per_threshold,
                                int min_thresh,
                                int max_thresh);

#endif // __ENRICHMENT_PROFILE__
