#ifndef __GENOME__
#define __GENOME__

/*
Encapsulation of a genome.

A genome is composed of multiple chromosomes. 

Each chromosome is added to a genome object.

Each chromosome can be queried for a specific nucleotide at a position, for its length, etc.
*/

#include <vector>
using namespace std;

// Following encoding ensures that the code for a nucleotide and its complement is 5.
#define NUC_COMPX (0) // Complement of X nucleotide.
#define NUC_A (1)
#define NUC_C (2)
#define NUC_G (3)
#define NUC_T (4)
#define NUC_X (5)

class t_chromosome;

// Initial form for the chromosome mappability map. 
class t_chromosome_mappability_map
{
public:
	t_chromosome_mappability_map(char* mapp_map_fp, t_chromosome* _chromosome); // Load the mappability map from a file.
	t_chromosome_mappability_map(t_chromosome* _chromosome); // Initialize the counts to 0 for all positions.
	~t_chromosome_mappability_map();

	t_chromosome* chromosome;
	int* multiplicity_counts; // These are the counts for how many times the corresponding nucleotide is repeated.
};

class t_chromosome
{
public:
	t_chromosome(char* chromosome_fp); // Initialize a genome object with one chromosome.
	~t_chromosome();

	char* nucs;
	int length;
	char* fp;

	// The quaternary codes for all of the ASCII characters that may be encountered in a fasta file.
	char* quaternary_codes_per_ASCII_char;
	void set_quaternary_codes();
	char complement(char nuc);

	// The nucleotide sequence for the genome.
	bool disk_mapped; // This flag sets up disk mapped data access, slow but does not load the whole genome into memory.

	// Access the nucleotide at certain chromosome at certain nucleotide position.
	char x_nuc(int i_nuc);
	char x_nuc_rev(int i_nuc);

	void load_chromosome(char* chromosome_fp);

	// Fill a fragment buffer (of coded nucleotides), using the chromosome information.
	bool fill_quaternary_coded_fragment_on_fore_strand(char* fragment, int l_fragment, int base_i);
	bool fill_quaternary_coded_fragment_on_rev_strand(char* fragment, int l_fragment, int base_i);

	// Go over the fragments and fill the buffer with the first valid one. Updates the base index automatically.
	bool fill_first_valid_quaternary_coded_fragment_on_fore_strand(char* fragment, int l_fragment, int& base_i);
	bool fill_first_valid_quaternary_coded_fragment_on_rev_strand(char* fragment, int l_fragment, int& base_i);

	// Converts (packs) four nucleotides to one byte. nucs is a pointer to (n_nucs) many nucleotides.
	static char pack_nucs_to_bin(char* nucs, int n_nucs);
	static void unpack_nucs_from_byte(char byte, char* nucs, int n_nucs);

	int l_chromosome();
};

#endif // __GENOME__


