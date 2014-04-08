#ifndef __CHROMOSOME_HASHTABLE__
#define __CHROMOSOME_HASHTABLE__

class t_mem_pool;
class t_chromosome;

// t_pos encapsulates a position in a chromosome, used for generating the hashtable for the chromosomes.
class t_pos
{
public:
	t_pos(int _index);
	~t_pos();
	int chr_index;
	t_pos* next;
	//bool found;
};

/*
Creates a simple hashtable with simply indexing fragments. No collisions are possible.

The hashtable is generated from one chromosome. It is basically a summarization of the information in the chromosome. 
This enables looking at another chromosome (or the same chromosome) to search and compare motifs.
*/
class t_chromosome_hashtable
{
public:
	t_chromosome_hashtable(char* _hashtable_fp); // Load the hashtable from file.
	t_chromosome_hashtable(t_chromosome* _chromosome, int _l_fragments);

	~t_chromosome_hashtable();

	// Memory pool to store the hash table. This is added by default, not #ifdef'ed any more.
	t_mem_pool* mem_pool;

	t_chromosome* chromosome; // Genome.
	int l_fragments; // Length of fragments that are indeced in hashtable.
	int n_fragments; // # of possible fragments.

	// Fragment size.
	void dump_hashtable(char* ht_dump_fp);
	void load_hashtable(char* _hashtable_fp);
	void set_hashtable();

	// These are the main hashtables. fore_strand_ht[index] is a pointer to a linked list of chromosome indices.
	// There are 4^(l_fragments) many t_pos*'s in this list. Each t_pos* is a pointer to the beginning of a linked list 
	// of chromosome indices which match the corresponding fragment.
	t_pos** fore_strand_ht;
	//t_pos** rev_strand_ht;

	// Following 4 functions are just for counting the number of possible fragments.
	bool init_fragment_counting_digits(char* frag_cnt);
	void update_fragment_counting_digits(char* frag_cnt);
	bool check_fragment_counting_end(char* frag_cnt);
	void buffer_fragment_counting_digit_string(char* buffer, char* frag_cnt);

	// Compute the hash for a fragment.
	int get_hash_by_quaternary_coded_fragment(char* fragment);
	int get_fore_str_hash_by_quaternary_coded_rev_str_fragment(char* fragment);

	// Compare the entries, this is utilized in debugging the hashtable dumping.
	bool compare(t_chromosome_hashtable* ht_2_compare);
	bool compare_recursive(t_pos* entry1, t_pos* entry2);

	// Find the fragment.
	void find_fragment(char* fragment, int n_max_mismatches);
	void get_consecutive_fragments(t_pos* l_frag_list, t_pos* r_frag_list);

	// Search a chromosome using this hashtable, mark the positions on the chromosome with how many times a certain fragment is repeated in the genome.
	// This process outputs a chromosome-wide mappability-map for the chromosome at hand. This is represented as a disk-mapped array.
	void create_mappability_map_per_chromosome(t_chromosome* chromosome_2_map, char* mappability_map_fp);

	void create_mappability_map_per_chromosome_10_nuc_fragment_ht(t_chromosome* chromosome_2_map, char* mappability_map_fp);
};

#endif // __CHROMOSOME_HASHTABLE__




