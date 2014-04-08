#ifndef __PEAKSEQ__
#define __PEAKSEQ__

//void peakseq(char* chr_list);
void peakseq(char* exp_id,
			char* chr_list_fp,
			int enrichment_fragment_length,
			char* background_model,
			double target_FDR,
			double max_Qvalue,
			int min_thresh,
			int max_thresh,
			int n_sims,
			int min_gap_per_bw_peaks,
			double P_f,
			vector<char*>* chip_seq_reads_data_dirs,
			vector<char*>* input_reads_data_dirs,
			char* mappability_map_fp,
			char* peak_prof_dump_dir,
			char* narrowPeak_op_fp,
			int simulation_seed);

#endif // __PEAKSEQ__
