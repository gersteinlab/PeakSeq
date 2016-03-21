#ifndef __INPUT_NORMALIZATION__
#define __INPUT_NORMALIZATION__

#include <vector>
using namespace std;

struct t_mapped_fragment;

double get_input_scaling_factor(int n_meg_wins,
                                vector<t_mapped_fragment*>* chip_seq_fragments,
                                vector<t_mapped_fragment*>* control_fragments,
                                double P_f,
                                int enrichment_mapped_fragment_length);

//double get_input_scaling_factor(t_chromosome* cur_chr, vector<t_mapped_fragment*>* cur_chr_chip_seq_fragments, vector<t_mapped_fragment*>* cur_chr_control_fragments, double P_f, int enrichment_mapped_fragment_length);
double pearson_correlation(vector<int>* control_frag_cnt, vector<int>* chip_seq_frag_cnt);
double PeakSeq_slope(vector<int>* control_frag_cnt, vector<int>* chip_frag_cnt, int n_processed_wins);
double slope(vector<int>* control_frag_cnt, vector<int>* chip_frag_cnt);

#endif // __INPUT_NORMALIZATION__
