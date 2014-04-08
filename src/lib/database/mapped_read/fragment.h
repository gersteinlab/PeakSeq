#ifndef __FRAGMENT__
#define __FRAGMENT__

#include <vector>
using namespace std;

const int MEG_BASE = 1000 * 1000;
const int K_BASE = 1000;

class t_chromosome;

#define FRAG_START (1)
#define FRAG_END (-1)

struct t_frag_end
{
	int pos;
	signed int side; // Is this an end or a beginning for the fragment? (Using the same trick as in current PeakSeq code.);
};

struct t_fragment
{
        int base_index;
        int sequenced_fragment_length; // Necessary when building enrichment profile. This is the fragment length coming from ChIP-Seq data.
        char strand_char;
};

// Fast version for fragment loading: Reads binary file, does not do validation of the fragment sequences.
void load_fragments_binary(char* mapped_reads_fp, vector<t_fragment*>* fore_strand_fragments, vector<t_fragment*>* rev_strand_frag);
void dump_fragments_binary(char* mapped_reads_bin_fp, vector<t_fragment*>* fragments);

void load_fragments(t_chromosome* chr, char* mapped_reads_fp, vector<t_fragment*>* fore_strand_fragments, vector<t_fragment*>* rev_strand_frag);
void delete_fragments(vector<t_fragment*>* fragment_list);
void delete_fragments(t_fragment** fragment_list);
bool check_fragment_quality(char* fragment, char* chr_subseq, char* quality_str);
char get_complement_nuc_by_nuc(char nuc_2_complement);
bool sort_fragments(t_fragment* frag1, t_fragment* frag2);
vector<t_fragment*>* prune_fragment_list_by_strand(vector<t_fragment*>* chr_fragments, char strand_2_prune, int n_max_reps);

/*
Following function is the only function that does tag extension with enrichment_fragment_length parameter. This should not be handled anywhere else.
*/
vector<t_fragment*>* forwardize_combine_sort_fore_rev_strand_frags(vector<t_fragment*>* fore_frag_list, vector<t_fragment*>* rev_frag_list, int enrichment_fragment_length);

int* get_n_reads_per_window(int n_wins, vector<t_fragment*>* frag_list);

vector<t_frag_end*>* get_frag_ends(vector<t_fragment*>* fragments);
void add_fragment_ends(vector<t_frag_end*>* frag_ends, int start, int length);
bool sort_ends(t_frag_end* end1, t_frag_end* end2);

void delete_frag_ends(vector<t_frag_end*>* frag_ends);

#endif // __FRAGMENT__
