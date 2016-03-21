#ifndef __OUTPUT__
#define __OUTPUT__

/*
Manages output to the file. (Singleton pattern)
*/
class t_peakseq_output
{
private:
	// Hidden constructor.
	t_peakseq_output(int i);

public:
	// Check if the object is allocated, if so do not do anything. Otherwise setup the data: Initialize the output file, copy the id.
	// Set the current chromosome id to emptry string.
	t_peakseq_output(char* _chip_seq_experiment_id);
	t_peakseq_output(char* _chip_seq_experiment_id, char* _narrowPeak_op_fp);
	~t_peakseq_output();

	static t_peakseq_output* singleton;
	char* cur_chr_id;
	char* chip_seq_experiment_id;

	static void set_cur_chr_id(char* new_chr_id);

	static void dump_PeakSeq_title(char* title_string);
	static void dump_PeakSeq_line(char* result_string); // Dump a result string.

	char* narrowPeak_op_fp;
	static void dump_narrowPeak_line(char* result_string);
};

#endif // __OUTPUT__
