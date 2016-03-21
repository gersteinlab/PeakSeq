#ifndef __PEAK_REGION__
#define __PEAK_REGION__

class t_peak_region
{
public:
	t_peak_region(int _start_i_nuc, int _end_i_nuc);
        t_peak_region(int _start_i_nuc, int _end_i_nuc, double log_p_val);

        t_peak_region();
	~t_peak_region();

	int start_i_nuc;
	int end_i_nuc;
	double log_p_val;
	double q_val;
	double enrichment;
	char* chr_id;

	double n_chip_seq_reads;
	double n_control_reads;
	double n_adj_control_reads;

	// Maximum height and index of maximum height in the peak.
	int max_height;
	int max_height_start_i;
	int max_height_end_i;
};

#endif // __PEAK_REGION__
