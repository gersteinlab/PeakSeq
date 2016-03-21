#ifndef __CHIP_SEQ_EXP_CHR_DATA__
#define __CHIP_SEQ_EXP_CHR_DATA__

#include <vector>
using namespace std;

struct t_mapped_fragment;
class t_enrichment_profile;

/*
This structure contains the information about a chromosome from a chip-seq experiment.
This is the data that is stored from the processed data files. Other parameters are read from
the command line.

Note that the chip seq data object does not include a reference to the chromosome. All the data validation 
should be done while writing the data files.
*/
class t_chip_seq_chr_data
{
public:
	// The id is used to form the file names associated with this chromosome.
	t_chip_seq_chr_data(vector<char*>* chip_seq_reads_data_dirs, 
						vector<char*>* input_reads_data_dirs, 
						char* _chr_id, 
						char* mappability_map_fp, 
						int enrichment_mapped_fragment_length);

	~t_chip_seq_chr_data();

    vector<t_mapped_fragment*>* cur_chr_chip_seq_fore_frags;
    vector<t_mapped_fragment*>* cur_chr_chip_seq_rev_frags;

	char* chr_id;

	// The set of chip seq and control fragments in the 
	vector<t_mapped_fragment*>* chip_seq_fragments;
	vector<t_mapped_fragment*>* control_fragments;

        // The enrichment profiles.
        t_enrichment_profile* chip_seq_enr_prof;
        t_enrichment_profile* control_enr_prof;

	// Number of uniquely mappable nucleotide counts per megabase window.
	vector<int>* n_uniquely_mappable_nucs_per_meg_win; 

	// Number of megabase windows in this chromosome: Use the 
	int n_meg_wins();

	// Load the mapability map file, which is optional.
	void load_mapability_map_file(char* mappability_map_fp);

	int n_nucleotides();
};

#endif // __CHIP_SEQ_EXP_DATA__


