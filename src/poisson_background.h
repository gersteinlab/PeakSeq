#ifndef __POISSON_BACKGROUND__
#define __POISSON_BACKGROUND__

#include <vector>
using namespace std;

class t_chromosome;
class t_enrichment_profile;
struct t_mapped_fragment;
class t_chip_seq_chr_data;

double* get_poisson_thresholds(t_chip_seq_chr_data* cur_chr_data,
								double target_fdr,
								int enrichment_mapped_fragment_length,
								int min_thresh,
								int max_thresh);

#endif // __POISSON_BACKGROUND__
