#ifndef __SIGNAL_TRACK_FILE_INTERFACE__
#define __SIGNAL_TRACK_FILE_INTERFACE__

#include <vector>

using namespace std;

/*
This the simple entry value format from the WIG file. This needs to be converted into sgr format before dumping into sgr format.
*/
struct t_wig_entry
{
	int pos;
	int span;
	double value;
};

struct t_profile_site
{
	int i_nuc;
	double height;
};

struct t_bedGraph_entry
{
	int start;
	int end;
	double value;
};

typedef vector<vector<t_profile_site*>*> t_signal_profile_list;

/*
t_signal_profile: Encapsulates a signal profile.
*/
struct t_signal_profile
{
	vector<char*>* chr_ids;
	t_signal_profile_list* profiles_per_chr;
};

t_signal_profile* allocate_signal_profile();

class t_peak_region;
struct t_annot_region;

bool sort_wig_entries(t_wig_entry* entry1, t_wig_entry* entry2);
bool sort_bedGraph_entries(t_bedGraph_entry* entry1, t_bedGraph_entry* entry2);
bool sort_profile_sites(t_profile_site* entry1, t_profile_site* entry2);

vector<t_annot_region*>* load_BED_with_avg_features(char* local_means_fp, int block_length);

void parse_WIG_info_line(char* wig_info_line, 
						char* chrom,
						int& start, 
						int& step, 
						int& span);

bool parse_new_WIG_entry(FILE* f_wig,
						vector<vector<t_wig_entry*>*>* sig_track_lists,
						vector<char*>* chr_ids,
						char* chrom, 
						int& start,
						int& step,
						int& span,
						int& signal_nuc_start, 
						double& signal_value);

bool parse_new_bedGraph_entry(FILE* f_bedGraph,
								vector<vector<t_wig_entry*>*>* sig_track_lists,
								vector<char*>* chr_ids,
								char* chrom, 
								int& start,
								int& step,
								int& span,
								int& signal_nuc_start, 
								double& signal_value);

// These functions are for preprocessing the signal track files: WIG, sgr, ... into the sgr format, which is the standard format used by t_enrichment_profile class.
// Load the WIG file and dump parsed SGR files. These can be loaded into enrichment profiles.
void parse_WIG_formatted_signal_track(vector<char*>* chr_fps, char* parsed_signal_tracks_op_dir, char* wig_fp);

void parse_bedGraph_formatted_signal_track(vector<char*>* chr_fps, char* parsed_signal_tracks_op_dir, char* wig_fp);

FILE* get_sig_track_file_pointer_by_chr_file_name(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_fn);
FILE* get_sig_track_file_pointer_by_chr_id(vector<FILE*>* fragment_f_ptrs, vector<char*>* chr_ids, char* chr_id);

vector<t_wig_entry*>* get_sig_track_list_pointer_by_chr_file_name(vector<vector<t_wig_entry>*>* sig_track_lists, vector<char*>* chr_ids, char* chr_fn);
vector<t_wig_entry*>* get_sig_track_list_pointer_by_chr_id(vector<vector<t_wig_entry*>*>* sig_track_lists, vector<char*>* chr_ids, char* chr_id);

void dump_sgr_per_WIG_entries(vector<t_wig_entry*>* wig_entries, FILE* f_sgr, char* chr_id);

vector<t_peak_region*>** get_peaks_per_thresholds(vector<t_profile_site*>* profile_sites, int min_threshold, int max_threshold, int min_gap_bw_peaks, double scaling_factor);

/*
Following are the t_profile_site functions, which provide the backend functionality for signal track processing.

Above functions should switch to t_profile_site versions. Basically what is needed is loading a wig file into t_profile_site list.

The chromosome handling is different from annot_regions. The reason for this is that t_profile_site entries do not contain the chromosome id's.
The loading of the files and the chromosome lists are done in parallel. This is similar to re-structuring of the annot_region lists in annot_region
interface. The reason for different chromosome handling is that the bedgraph files usually contain very large number of entries and it is very 
costly to get sublist of profiles for each chromosome every time the list is processed. Therefore, the chromosome list is stored with the separate signal
profiles.
*/
void dump_peak_profiles_per_list(vector<t_profile_site*>* profile_sites, 
									vector<t_peak_region*>* peak_list, 
									char* dump_dir, 
									char* file_name_fmt_str);

void dump_profile(vector<t_profile_site*>* profile_sites,
					char* fp,
					char* chr_id,
					int min_i_nuc,
					int max_i_nuc, 
					double scaling_factor);

void dump_profile(t_signal_profile* signal_profile,
					char* fp,
					int min_i_nuc, 
					int max_i_nuc, 
					double scaling_factor);

void dump_profile_per_list(vector<t_profile_site*>* profile_sites,
							char* dump_dir,
							vector<t_annot_region*>* annot_regions,
							double scaling_factor);

void reverse_profile(vector<t_profile_site*>* profile_sites);
vector<t_profile_site*>* subtract_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2);
void square_profile(vector<t_profile_site*>* profile);
vector<t_profile_site*>* aggregate_profiles(vector<vector<t_profile_site*>*>* profile_lists);
vector<t_profile_site*>* get_variance_profile(vector<vector<t_profile_site*>*>* profile_list);

// Assumption is that the profile sites and annot_regions are in the same chromosome. The chromosome handling must be donw beforehand.
t_signal_profile_list* buffer_profile_per_list(vector<t_profile_site*>* profile_sites,
													vector<t_annot_region*>* annot_regions,
													double scaling_factor);

vector<t_profile_site*>* buffer_profile_per_region(vector<t_profile_site*>* profile_sites,
													int start, int end,
													double scaling_factor);

vector<t_profile_site*>* get_profile_per_delta_profile(vector<t_profile_site*>* delta_profile);
vector<t_profile_site*>* get_delta_profile_per_profile(vector<t_profile_site*>* profile);

vector<t_profile_site*>* aggregate_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2);

void translate_profile(vector<t_profile_site*>* profile, int new_start);

void amplify_profile(vector<t_profile_site*>* profile, double amp);

// Compute the simple pearson correlation between two tracks.
double correlate_profiles(vector<t_profile_site*>* profile1, vector<t_profile_site*>* profile2, int bin_size);

// Correlate genome wide profiles by taking the chromosomes into consideration, differs only in the bin generation.
double correlate_profiles(t_signal_profile* profile1, t_signal_profile* profile2, int bin_size);

double integrate_signal(vector<t_profile_site*>* profile);

// Dump the average signalper windows for profile sites: Signal aggregation over the windows.
//void dump_avg_signal_per_windows(vector<char*>* chr_ids, vector<vector<t_profile_site*>*>* profile_sites, char* fp, int win_size);
void dump_avg_signal_per_windows(vector<char*>* chr_ids, 
								t_signal_profile* signal_profile,
								char* fp, 
								int win_size);

void dump_binary_profile(t_signal_profile* signal_profile, 
						char* op_fp);

// Open a bedgraph file and load only the chromosome ids.
vector<char*>* get_bedGraph_chr_ids(char* bedGraph_fp);
vector<char*>* get_sgr_chr_ids(char* bedGraph_fp);

//void load_sgr(char* sgr_fp, 
//				vector<vector<t_profile_site*>*>* profile_sites_per_chroms,
//				vector<char*>* chr_ids);

//void load_bedGraph(char* bgr_fp, 
//				vector<vector<t_profile_site*>*>* profile_sites_per_chrom, 
//				vector<char*>* chrom);

t_signal_profile* load_sgr(char* sgr_fp);

t_signal_profile* load_bedGraph(char* bgr_fp);

void normalize_signal_profile_length(vector<t_profile_site*>* sig_prof, int total_normalized_width);

vector<double>* get_integrated_signals_per_bins(vector<t_profile_site*>* sites, t_annot_region* annot_region, int bin_size);

void delete_profile_sites(vector<t_profile_site*>* profile_sites);

void delete_signal_profile(t_signal_profile* signal_profile);

// Translate all the profiles to leftmost position.
void translate_profiles_to_origin(vector<vector<t_profile_site*>*>* profiles);

double correlate(vector<double>* sig1, vector<double>* sig2);
double correlate_non_zeros(vector<double>* bgr1_sigs_per_bin, vector<double>* bgr2_sigs_per_bin);

// Binary file loading interfaces.
t_signal_profile* load_binary_profile(char* fp);
void binarize_bedGraph(char* fp, char* op_fp);

#endif // __SIGNAL_TRACK_FILE_INTERFACE__

