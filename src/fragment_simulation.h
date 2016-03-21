#ifndef __FRAGMENT_SIMULATION__
#define __FRAGMENT_SIMULATION__

//void generate_simulated_fragments(char* bed_fp, int n_frags_per_region, int l_frag, char* chip_op_tagalign_fp, char* ip_op_tagalign_fp);
void generate_simulated_fragments(char* bed_fp, int n_chip_frags_per_region, int n_ip_frags_per_region, int l_frag, char* chip_op_tagalign_fp, char* ip_op_tagalign_fp);
void get_uniform_cdf_per_region(int l_region, vector<double>* cdf);

#endif // __FRAGMENT_SIMULATION__