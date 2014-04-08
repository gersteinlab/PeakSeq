#include <stdio.h>
#include <stdlib.h>
#include "chromosome_hashtable.h"
#include "chromosome.h"
#include "../utils/memory/mem_pool.h"
#include "../utils/file/utils.h"
#include <string.h>

t_pos::t_pos(int _chr_index)
{
	this->chr_index = _chr_index;
	this->next = NULL;
}

t_pos::~t_pos()
{
	if(this->next != NULL)
	{
		delete(this->next);
	}
}

// Construct an empty hashtable to be loaded later. Otherwise it is not of any use.
t_chromosome_hashtable::t_chromosome_hashtable(char* _hashtable_fp)
{
	this->chromosome = NULL;
	this->l_fragments = 0; // length of the fragments that are used to index the table. The size of the table is 4^(l_fragments).
	this->n_fragments = 0;
	this->mem_pool = NULL;

	// Load from the file. Note that this hashtable object will not have genome and chromosome information, which is not really needed 
	// because, in theory, the genome (corresponding to the hashtable) can be reconstructed using the hashtable alone. Adding genome here 
	// creates redundancy.
	this->load_hashtable(_hashtable_fp);
}

// Genome can be NULL, in the case that hashtable is going to be loaded.
t_chromosome_hashtable::t_chromosome_hashtable(t_chromosome* _chromosome, int _l_fragments)
{
	this->chromosome = _chromosome;
	this->l_fragments = _l_fragments; // length of the fragments that are used to index the table. The size of the table is 4^(l_fragments).
	this->mem_pool = NULL;

	// Count the fragments.
	char* frag_digits = new char[this->l_fragments + 3];
	//this->init_fragment_counting_digits(frag_digits);
	//this->n_fragments = 1; // Set to 1 to count the inited value in previous expression.
	//do
	//{
	//	this->update_fragment_counting_digits(frag_digits);
	//	this->n_fragments++;
	//	//printf("%d             \r", this->n_fragments);

	//	//if(this->n_fragments - 1 != this->get_hash_by_quaternary_coded_fragment(frag_digits))
	//	//{
	//	//	char buffer[30];
	//	//	this->buffer_fragment_counting_digit_string(buffer, frag_digits);
	//	//	printf("Counting failed: %d, %d (%s)\n", this->n_fragments - 1, this->get_hash_by_quaternary_coded_fragment(frag_digits), buffer);
	//	//	exit(0);
	//	//}
	//}
	//while(!this->check_fragment_counting_end(frag_digits));

	char* max_fragment = new char[this->l_fragments];
	for(int i_nuc = 0; i_nuc < this->l_fragments; i_nuc++)
	{
		max_fragment[i_nuc] = NUC_T;
	} // i_nuc loop.

	this->n_fragments = this->get_hash_by_quaternary_coded_fragment(max_fragment) + 2;

	printf("Counted %d fragments with %d nucleotides.\n", this->n_fragments, this->l_fragments);

	delete[] frag_digits;

	this->set_hashtable();
}

t_chromosome_hashtable::~t_chromosome_hashtable()
{
	// Note that mem_pool includes all the memory for the hastable structure.
	if(this->mem_pool != NULL)
	{
		delete(this->mem_pool);
	}
	else
	{
		// Must manually free the hashtable.
		for(int i_frag = 0; i_frag < this->n_fragments; i_frag++)
		{
			// Following calls recursively deletion of all the nodes in the linked list.
			delete(this->fore_strand_ht[i_frag]);
		} // i_frag loop.
	}

	// ht array is allocated from heap, not from mem_pool.
	delete [] this->fore_strand_ht;
}

void t_chromosome_hashtable::dump_hashtable(char* ht_dump_fp)
{
	//char ht_dump_fp[1000];
	//sprintf(ht_dump_fp, "%s.ht", this->genome->chromosomes->at(this->i_chromosome)->fp);
	printf("Dumping hashtable to %s.\n", ht_dump_fp);

	FILE* f_ht = fopen(ht_dump_fp, "wb");
	if(f_ht == NULL)
	{
		printf("Could not open %s for writing.\n", ht_dump_fp);
		exit(0);
	}

	// Add the length of fragments in the hashtable.
	fwrite(&(this->l_fragments), sizeof(int), 1, f_ht);

	for(int cur_hash = 0; cur_hash < this->n_fragments; cur_hash++)
	{
		// Dump the entries for the current fragment entry. 
		t_pos* cur_pos_entry = this->fore_strand_ht[cur_hash];
		while(cur_pos_entry != NULL)
		{
			fwrite(&(cur_pos_entry->chr_index), sizeof(int), 1, f_ht);
			cur_pos_entry = cur_pos_entry->next;
		} // go over all the t_pos elements for this hashtable entry.

		// Mark the end of this entry.
		int i_end = 0xffffffff;
		fwrite(&(i_end), sizeof(int), 1, f_ht);
	} // cur_hash loop.

	fclose(f_ht);
}


void t_chromosome_hashtable::set_hashtable()
{
	if(this->chromosome == NULL)
	{
		printf("Cannot create a hastable without chromosome data.\n");
		return;
	}

	printf("Building the hashtable.\n");

	this->fore_strand_ht = new t_pos*[this->n_fragments+3];

	for(int i = 0; i < this->n_fragments; i++)
	{
		this->fore_strand_ht[i] = NULL; // Initialize to NULL.
	} // i loop

	// Populate the hashtable: 
	char* cur_fragment = new char[this->l_fragments + 3];

	// Initial fill for fragment: The first index to the nucleotides in the chromosome is 1.
	int base_i = 1;

	// Initialize the memory pool.
	this->mem_pool = new t_mem_pool();

	// Fill the first valid fragment on the forward strand.
	while(this->chromosome->fill_first_valid_quaternary_coded_fragment_on_fore_strand(cur_fragment, this->l_fragments, base_i))
	{
		//printf("Base_i: %s                         \r", base_i, cur_fragment);
		//getc(stdin);

/*
        if(base_i % 10000 == 0)
        {
			FILE* f_dump = fopen("dump.txt", "a");
			fprintf(f_dump, "%d(%d)\n", base_i, this->chromosome->l_chromosome() - this->l_fragments);
			fclose(f_dump);
        }
*/
//		char buff[15];
//		this->buffer_fragment_counting_digit_string(buff, cur_fragment);
//		printf("%s    \n", buff);

		// Get the hash for this fragment and add this position.
		int fragment_index = this->get_hash_by_quaternary_coded_fragment(cur_fragment);

		//printf("Hash is %d\n", fragment_index);
		//getc(stdin);

		// Add a new node to the beginning of the list.
		// Store the node that was previously at the beginning.
		t_pos* cur_leader_node = this->fore_strand_ht[fragment_index];

		// Allocate a new beginning node. The nodes are sorted from 3'-most to 5'-most.
		this->fore_strand_ht[fragment_index] = (t_pos*)this->mem_pool->get_mem_from_pool(sizeof(t_pos));
		this->fore_strand_ht[fragment_index]->chr_index = base_i;
		this->fore_strand_ht[fragment_index]->next = cur_leader_node; // Set the next to the previous head node.

		// Move to the next fragment.
		base_i++;
	} // i_nuc loop

	delete[] cur_fragment;
}

/*
Following REVERSES the ordering of the fragments in the linked lists.
*/
void t_chromosome_hashtable::load_hashtable(char* _hashtable_fp)
{
	printf("Loading hashtable from %s.\n", _hashtable_fp);

	FILE* f_ht = open_f(_hashtable_fp, "rb");

	// Read the length of fragments.
	fread(&this->l_fragments, sizeof(int), 1, f_ht);

	printf("l_fragments is read as %d\n", this->l_fragments);

	// Count the number of fragments (the usual way).
	//char* frag_digits = new char[this->l_fragments + 3];
	//this->init_fragment_counting_digits(frag_digits);
	//this->n_fragments = 1;
	//do
	//{
	//	this->update_fragment_counting_digits(frag_digits);
	//	this->n_fragments++;
	//}
	//while(!this->check_fragment_counting_end(frag_digits));

	char* max_fragment = new char[this->l_fragments];
	for(int i_nuc = 0; i_nuc < this->l_fragments; i_nuc++)
	{
		max_fragment[i_nuc] = NUC_T;
	} // i_nuc loop.

	this->n_fragments = this->get_hash_by_quaternary_coded_fragment(max_fragment) + 2;

	printf("Counted %d fragments with %d nucleotides.\n", this->n_fragments, this->l_fragments);

	//delete [] frag_digits;
	delete [] max_fragment;

	// It is necessary to compute n_fragments before doing following. 
        this->mem_pool = new t_mem_pool();

        // Allocate and populate the second hashtable.
		this->fore_strand_ht = new t_pos*[this->n_fragments + 2];

        for(int i_frag = 0; i_frag < this->n_fragments; i_frag++)
        {
                this->fore_strand_ht[i_frag] = NULL;
        } // i_frag loop.


	// Following loading reverses the order in which the fragments were dumped.
	int i_frag = 0;
	while(1)
	{
		int cur_val = 0;
		if(fread(&cur_val, sizeof(int), 1, f_ht) != 1)
		{
			break;
		}

		if(cur_val == 0xffffffff)
		{
			// Move to the next entry in the hashtable.
			i_frag++;
		}
		else
		{
				t_pos* leader_node = this->fore_strand_ht[i_frag];

				this->fore_strand_ht[i_frag] = (t_pos*)this->mem_pool->get_mem_from_pool(sizeof(t_pos));
				this->fore_strand_ht[i_frag]->chr_index = cur_val;
				this->fore_strand_ht[i_frag]->next = leader_node;
		}
	} // hash table file read loop.

	fclose(f_ht);

	// Verify the monotonous INCREASING chr_index for nodes in the same list.
	// The loaded list is monotonously increasing because the loading loop above 
	// reverses the monotonously decreasing order of nodes in the saved list.
	for(i_frag = 0; i_frag < this->n_fragments; i_frag++)
    {
		t_pos* cur_node = this->fore_strand_ht[i_frag];
		t_pos* cur_next_node = NULL;	
		while(cur_node != NULL && 
				cur_node->next != NULL)
		{
			cur_next_node = cur_node->next;

			// Is cur_node placed after cur_next_node in the genome? This should not happen.
			if(cur_node->chr_index > cur_next_node->chr_index)
			{
				printf("Problematic positioning of nodes in the loaded table.\n");
				exit(0);
			}

			cur_node = cur_node->next;
		}
	} // i_frag loop.

}

// Compute the hash for a fragment.
int t_chromosome_hashtable::get_hash_by_quaternary_coded_fragment(char* fragment)
{
	int hash = 0;
	int base = 1; // 4^0.

	// 5' nucleotide is the least significant nucleotide.
	for(int i_nuc = 0; i_nuc < this->l_fragments; i_nuc++)
	{
		hash += base * (fragment[i_nuc] - 1);
		base *= 4;
	}

	return(hash);
}

/*
Following function enables utilizing the forward strand hashtable for reverse strand: Determine the forward strand
hash from the reverse strand fragment. This is useful because when searching reverse strand, the fragment that is being
searched needs to be mapped to forward strand to find the complementary fragments.
*/
int t_chromosome_hashtable::get_fore_str_hash_by_quaternary_coded_rev_str_fragment(char* fragment)
{
	int hash = 0;
	int base = 1; // 4^0.

	// Start from the 3' side, complement the nucleotides and sum up.
	for(int i_nuc = this->l_fragments - 1; i_nuc >= 0; i_nuc--)
	{
		// Following uses the shortcut to complement nucleotides.
		hash += base * ((5 - fragment[i_nuc]) - 1);
		base *= 4;
	}

	return(hash);
}


// Compare the entries, this is utilized in debugging the hashtable dumping.
bool t_chromosome_hashtable::compare_recursive(t_pos* entry1, t_pos* entry2)
{
	if(entry1 == NULL)
	{
		if(entry2 == NULL)
		{
			return(true);
		}
		else
		{
			return(false);
		}
	}
	else if(entry1->chr_index != entry2->chr_index)
	{
		return(false);
	}
	else
	{
		if(entry1->next == NULL)
		{
			if(entry2->next == NULL)
			{
				return(true);
			}
			else
			{
				return(false);
			}
		}
		else
		{
			if(entry2->next == NULL)
			{
				return(false);
			}
			else
			{
				return(compare_recursive(entry1->next, entry2->next));
			}
		}
	}
}

bool t_chromosome_hashtable::init_fragment_counting_digits(char* frag_cnt)
{
	for(int i = 0; i < this->l_fragments; i++)
	{
		frag_cnt[i] = NUC_A;
	} // i loop

	return(true);
}

void t_chromosome_hashtable::update_fragment_counting_digits(char* frag_cnt)
{
	// Go over all nucleotides in the fragment and update them starting from the least significant one.
	bool carry = false;
	for(int i = this->l_fragments-1; i >= 0; i--)
	//for(int i = 0; i < this->l_fragments; i++)
	{	
		// Increment the least significant digit by one.
		if(i == this->l_fragments-1)
		{
			// Include a carry if necessary.
			if(frag_cnt[i] == NUC_T)
			{
				frag_cnt[i] = NUC_A;
				carry = true;
			}
			else
			{
				frag_cnt[i]++;
			}
		}
		else
		{
			// Process carry.
			if(carry)
			{
				if(frag_cnt[i] == NUC_T)
				{
					frag_cnt[i] = NUC_A;
					carry = true; // Push the carry upward.
				}
				else
				{
					// Carry processed, no further carries are necessary.
					frag_cnt[i]++;
					carry = false;
				}
			}
			else
			{
				// There is no carry to process. Nothing has to be done for upward digits.
			} // carry check.
		}
	} // i loop
}

bool t_chromosome_hashtable::check_fragment_counting_end(char* frag_cnt)
{
	for(int i = 0; i < this->l_fragments; i++)
	{
		if(frag_cnt[i] != NUC_T)
		{
			return(false);
		}
	} // i loop

	return(true);
}

void t_chromosome_hashtable::buffer_fragment_counting_digit_string(char* buffer, char* frag_cnt)
{
	char nucs[] = "NACGTN";
	for(int i = 0; i < this->l_fragments; i++)
	{
		//printf("%c", nucs[frag_cnt[i]]);
		buffer[i] = nucs[frag_cnt[i]];
	} // i loop

	buffer[l_fragments] = 0;
}

bool t_chromosome_hashtable::compare(t_chromosome_hashtable* ht_2_compare)
{
	// Compare the hashtables.
	for(int i_frag = 0; i_frag < this->n_fragments; i_frag++)
	{
		printf("Comparing %d. entry.\r", i_frag);

		// Compare the position entries.
		t_pos* ht_pos_entry = this->fore_strand_ht[i_frag];

		t_pos* ht_2_pos_entry = ht_2_compare->fore_strand_ht[i_frag];

		if(!compare_recursive(ht_pos_entry, ht_2_pos_entry))
		{
			return(false);
		}
	} // i_frag loop for htable comparison.	
}

/*
Find a fragment with maximum number of mismatches.
*/
void t_chromosome_hashtable::find_fragment(char* fragment, int n_max_mismatches)
{
	if(strlen(fragment) != (2 * this->l_fragments))
	{
		printf("Cannot search for a fragment with %d nucleotides in a hashtable of fragments of length %d.\n", strlen(fragment), this->l_fragments);
		exit(0);
	}	

	char* l_frag = new char[this->l_fragments + 3];
	char* r_frag = new char[this->l_fragments + 3];

	// Need to re-encode the fragments in quaternary coding.
	for(int i = 0; i < this->l_fragments; i++)
	{
		l_frag[i] = this->chromosome->quaternary_codes_per_ASCII_char[fragment[i]];
	}
	for(int i = this->l_fragments; i < strlen(fragment); i++)
	{
		r_frag[i - this->l_fragments] = this->chromosome->quaternary_codes_per_ASCII_char[fragment[i]];
	}

	printf("Searching forward strand for %s\n", fragment, l_frag, r_frag);

	// Search for the forward strand.
	int l_frag_hash = this->get_hash_by_quaternary_coded_fragment(l_frag);
	int r_frag_hash = this->get_hash_by_quaternary_coded_fragment(r_frag);

	// Get the positions where these are placed in the other chromosome.
	t_pos* l_frag_list = this->fore_strand_ht[l_frag_hash];
	t_pos* r_frag_list = this->fore_strand_ht[r_frag_hash];

	this->get_consecutive_fragments(l_frag_list, r_frag_list);

	printf("Searching backward strand for %s\n", fragment);

	// Treat the left and right fragment as if they are on reverse strand. 
	// Note that the order is changed: l_frag_hash is determined from right fragment and
	// r_frag_hash is determined from the left fragment.
	l_frag_hash = this->get_fore_str_hash_by_quaternary_coded_rev_str_fragment(r_frag);
	r_frag_hash = this->get_fore_str_hash_by_quaternary_coded_rev_str_fragment(l_frag);
	
	// The remaining search is the same as forward strand search.

	// Get the positions where these are placed in the other chromosome.
	l_frag_list = this->fore_strand_ht[l_frag_hash];
	r_frag_list = this->fore_strand_ht[r_frag_hash];

	this->get_consecutive_fragments(l_frag_list, r_frag_list);

	delete [] l_frag;
	delete [] r_frag;
}

void t_chromosome_hashtable::get_consecutive_fragments(t_pos* l_frag_list, t_pos* r_frag_list)
{
	t_pos* next_r_pos_entry = r_frag_list;

	// Trace all the nodes in order.
	t_pos* cur_l_entry = l_frag_list;

	// Trace all the nodes in order.
	cur_l_entry = l_frag_list;
	while(cur_l_entry != NULL && 
		next_r_pos_entry != NULL)
	{
		// Find the first right node that has the index higher than current_l_entry.
		while(next_r_pos_entry != NULL &&
				next_r_pos_entry->chr_index < cur_l_entry->chr_index)
		{
			next_r_pos_entry = next_r_pos_entry->next;
		}

		// Check to the rightmost vicinity of cur_l_node to see if there is a consecutive right node.
		t_pos* cur_next_r_node = next_r_pos_entry;
		bool search_done = false;
		while(!search_done &&
			cur_next_r_node != NULL)
		{
			if(cur_next_r_node->chr_index == cur_l_entry->chr_index + this->l_fragments)
			{
				search_done = true;
				printf("Found a copy @ %d\n", cur_l_entry->chr_index);
			}
			else if(cur_next_r_node->chr_index > cur_l_entry->chr_index + this->l_fragments)
			{
				// Right node just got out of range of search, it is not necessary to search any more.
				search_done = true;
			}
			cur_next_r_node = cur_next_r_node->next;
		}

		// Move to the next left entry.
		cur_l_entry = cur_l_entry->next;
	}
}

/*
Note that in order to run this function, it is only necessary to load a genome (disk mapped) and a
chromosome hashtable. Note that there is no expectation that the mapping is going for chromosomes in the same genome.
*/
void t_chromosome_hashtable::create_mappability_map_per_chromosome(t_chromosome* chromosome_2_map, char* mappability_map_fp)
{
	FILE* f_mappability_map = open_f(mappability_map_fp, "wb");

	printf("Creating the mappability map for chromosome %s into %s.\n", chromosome_2_map, mappability_map_fp);
	int base_i = 0;

	// Move over the chromosome_2_map.
	char cur_fragment[35];

	// READ FROM THE GENOME TO MAP and SEARCH OVER THE HASHTABLE, WHICH IS A SEARCH-EFFICIENT REPRESENTATION OF THE ORIGINAL GENOME.
	while(chromosome_2_map->fill_first_valid_quaternary_coded_fragment_on_fore_strand(cur_fragment, 30, base_i))
	{
		//printf("base_i = %d ", base_i);
		//fflush(stdout);
		char l_frag[30];
		char r_frag[30];

		//printf("Copying fragments.\n");
		//fflush(stdout);

		if(this->l_fragments != 15)
		{
			printf("The mappability map cannot be created with a hashtable of fragment size %d, need 15 long fragments.\n", this->l_fragments);
			exit(0);
		}

		// THE FRAGMENT LENGTHS ARE HARDCODED INTO HERE, THIS SHOULD BE CHANGED!!!!
		for(int i = 0; i < 15; i++)
		{
			l_frag[i] = cur_fragment[i];
		}
		for(int i = 15; i < 30; i++)
		{
			r_frag[i - 15] = cur_fragment[i];
		}

		int n_cur_frag_copies = 0;

/*
		if(base_i % 10000 == 0)
		{
			FILE* f_dump = fopen("dump.txt", "a");
			fprintf(f_dump, "%d\n", base_i);
			fclose(f_dump);
		}
*/
		//char buff[30];
		//this->buffer_fragment_counting_digit_string(buff, l_frag);
		//printf("l_frag: %s\n", buff);
		//this->buffer_fragment_counting_digit_string(buff, m_frag);
		//printf("m_frag: %s\n", buff);
		//this->buffer_fragment_counting_digit_string(buff, r_frag);
		//printf("r_frag: %s\n", buff);

		// Compute the hashes for all fragments from the chromosome to be mapped.
		int l_frag_hash = this->get_hash_by_quaternary_coded_fragment(l_frag);
		//int m_frag_hash = this->get_hash_by_quaternary_coded_fragment(m_frag);
		int r_frag_hash = this->get_hash_by_quaternary_coded_fragment(r_frag);

		// Get the positions where these are placed in the other chromosome.
		t_pos* l_frag_list = this->fore_strand_ht[l_frag_hash];
		//t_pos* m_frag_list = this->fore_strand_ht[m_frag_hash];
		t_pos* r_frag_list = this->fore_strand_ht[r_frag_hash];

		//printf("Computing sizes.\n");
		//fflush(stdout);
/*
		t_pos* l_start_node_ptr = l_frag_list;
		int l_size = 0;	
		t_pos* l_node = l_frag_list;
                int min_l_chr_index = l_node->chr_index;
                int max_l_chr_index = l_node->chr_index;
		while(l_node != NULL)
		{
			l_node = l_node->next;
			l_size++;
		}	
		printf("(%d, ", l_size);
*/
/*
                int m_size = 0;
                t_pos* m_node = m_frag_list;
                while(m_node != NULL)
                {
                        m_node = m_node->next;
                        m_size++;
                }
                printf("%d, ", m_size);
*/

/*
                int r_size = 0;
                t_pos* r_node = r_frag_list;
                int min_r_chr_index = r_node->chr_index;
                int max_r_chr_index = r_node->chr_index;
                while(r_node != NULL)
                {
                        r_node = r_node->next;
                        r_size++;
                }
                printf("%d) ", r_size);

		if(l_size > 1000 || r_size > 1000)
		{
			char buffer[40];
			buffer_fragment_counting_digit_string(buffer, cur_fragment);
			printf("%d: %d, %d (%s)\n", base_i, l_size, r_size, buffer);
			getc(stdin);
		}
*/
		// Go over all the nodes.
		// This is the node that is just right to the current left node.
		t_pos* next_r_pos_entry = r_frag_list;

		// Trace all the nodes in order.
		t_pos* cur_l_entry = l_frag_list;
		while(cur_l_entry != NULL && 
			next_r_pos_entry != NULL)
		{
			// Find the first right node that has the index higher than current_l_entry.
			while(next_r_pos_entry != NULL &&
					next_r_pos_entry->chr_index < cur_l_entry->chr_index)
			{
				next_r_pos_entry = next_r_pos_entry->next;
			}

			// Check to the rightmost vicinity of cur_l_node to see if there is a consecutive right node.
			t_pos* cur_next_r_node = next_r_pos_entry;
			bool search_done = false;
			while(!search_done &&
				cur_next_r_node != NULL)
			{
				if(cur_next_r_node->chr_index == cur_l_entry->chr_index + this->l_fragments)
				{
					search_done = true;
					n_cur_frag_copies++;
				}
				else if(cur_next_r_node->chr_index > cur_l_entry->chr_index + this->l_fragments)
				{
					// Right node just got out of range of search, it is not necessary to search any more.
					search_done = true;
				}
				cur_next_r_node = cur_next_r_node->next;
			}

			// Move to the next left entry.
			cur_l_entry = cur_l_entry->next;
		}
	
		//printf("[%d]         \n", n_cur_frag_copies);
	
		// Dump the count for this position in the mappability file. Note that most probably this will not be needed
		// later, only thing that is necessary is whether this fragment is uniquely mappable or not.
		fwrite(&n_cur_frag_copies, sizeof(int), 1, f_mappability_map);

		base_i++; // Update the position.
    }  //base_i loop.

	// Do the counts on the reverse strand of the chromosome for which 

	fclose(f_mappability_map);
}


void t_chromosome_hashtable::create_mappability_map_per_chromosome_10_nuc_fragment_ht(t_chromosome* chromosome_2_map, char* mappability_map_fp)
{
	FILE* f_mappability_map = open_f(mappability_map_fp, "wb");

	printf("Creating the mappability map for chromosome %s into %s.\n", chromosome_2_map, mappability_map_fp);
	int base_i = 0;

	// Move over the chromosome_2_map.
	char cur_fragment[35];

	if(this->l_fragments != 10)
	{
		printf("Hash table fragment length is not 10 nucleotides, cannot compute mappability with create_mappability_map_per_chromosome_10_nuc_fragment_ht when fragment length is not 10 nucleotides.\n");
		exit(0);
	}

	// READ FROM THE GENOME TO MAP and SEARCH OVER THE HASHTABLE, WHICH IS A SEARCH-EFFICIENT REPRESENTATION OF THE ORIGINAL GENOME.
	while(chromosome_2_map->fill_first_valid_quaternary_coded_fragment_on_fore_strand(cur_fragment, 30, base_i))
	{
		printf("base_i = %d ", base_i);
		fflush(stdout);

		char l_frag[30];
		char m_frag[30];
		char r_frag[30];

		for(int i = 0; i < 10; i++)
		{
			l_frag[i] = cur_fragment[i];
		}
		for(int i = 10; i < 20; i++)
		{
			m_frag[i - 10] = cur_fragment[i];
		}
		for(int i = 20; i < 30; i++)
		{
			r_frag[i - 20] = cur_fragment[i];
		}

		int n_cur_frag_copies = 0;

		//printf("Copying fragments.\n");
		//fflush(stdout);

		//if(base_i % 10000 == 0)
		//{
		//	FILE* f_dump = open_f("dump.txt", "a");
		//	fprintf(f_dump, "%d\n", base_i);
		//	fclose(f_dump);
		//}

		//char buff[30];
		//this->buffer_fragment_counting_digit_string(buff, l_frag);
		//printf("l_frag: %s\n", buff);
		//this->buffer_fragment_counting_digit_string(buff, m_frag);
		//printf("m_frag: %s\n", buff);
		//this->buffer_fragment_counting_digit_string(buff, r_frag);
		//printf("r_frag: %s\n", buff);

		// Compute the hashes for all fragments from the chromosome to be mapped.
		int l_frag_hash = this->get_hash_by_quaternary_coded_fragment(l_frag);
		int m_frag_hash = this->get_hash_by_quaternary_coded_fragment(m_frag);
		int r_frag_hash = this->get_hash_by_quaternary_coded_fragment(r_frag);

		// Get the positions where these are placed in the other chromosome.
		t_pos* l_frag_list = this->fore_strand_ht[l_frag_hash];
		t_pos* m_frag_list = this->fore_strand_ht[m_frag_hash];
		t_pos* r_frag_list = this->fore_strand_ht[r_frag_hash];

		//printf("Computing sizes.\n");
		//fflush(stdout);
/*
		t_pos* l_start_node_ptr = l_frag_list;
		int l_size = 0;	
		t_pos* l_node = l_frag_list;
                int min_l_chr_index = l_node->chr_index;
                int max_l_chr_index = l_node->chr_index;
		while(l_node != NULL)
		{
			l_node = l_node->next;
			l_size++;
		}	
		printf("(%d, ", l_size);
*/
/*
                int m_size = 0;
                t_pos* m_node = m_frag_list;
                while(m_node != NULL)
                {
                        m_node = m_node->next;
                        m_size++;
                }
                printf("%d, ", m_size);
*/

/*
                int r_size = 0;
                t_pos* r_node = r_frag_list;
                int min_r_chr_index = r_node->chr_index;
                int max_r_chr_index = r_node->chr_index;
                while(r_node != NULL)
                {
                        r_node = r_node->next;
                        r_size++;
                }
                printf("%d) ", r_size);

		if(l_size > 1000 || r_size > 1000)
		{
			char buffer[40];
			buffer_fragment_counting_digit_string(buffer, cur_fragment);
			printf("%d: %d, %d (%s)\n", base_i, l_size, r_size, buffer);
			getc(stdin);
		}
*/
		// Go over all the nodes.
		// This is the node that is just right to the current left node.
		t_pos* next_r_pos_entry = r_frag_list;
		t_pos* next_m_pos_entry = m_frag_list;

		// Trace all the nodes in order.
		t_pos* cur_l_entry = l_frag_list;
		while(cur_l_entry != NULL && 
			next_m_pos_entry != NULL && 
			next_r_pos_entry != NULL)
		{
			// Find the first middle node to the right of current left node.
			while(next_m_pos_entry != NULL &&
					next_m_pos_entry->chr_index < cur_l_entry->chr_index)
			{
				next_m_pos_entry = next_m_pos_entry->next;
			}

			// Find the first right node to the right of current middle node.
			while(next_m_pos_entry != NULL &&
					next_r_pos_entry != NULL &&
					next_r_pos_entry->chr_index < next_m_pos_entry->chr_index)
			{
				next_r_pos_entry = next_r_pos_entry->next;
			}

			// Check to the rightmost vicinity of cur_l_node to see if there is a consecutive middle node.
			t_pos* cur_next_m_node = next_m_pos_entry;
			bool search_done = false;
			while(!search_done &&
				cur_next_m_node != NULL)
			{
				if(cur_next_m_node->chr_index == cur_l_entry->chr_index + this->l_fragments)
				{
					// Check to the rightmost vicinity of cur_next_m_node to see if there is a consecutive right node.
					t_pos* cur_next_r_node = next_r_pos_entry;
					while(!search_done &&
						cur_next_r_node != NULL)
					{
						if(cur_next_r_node->chr_index == cur_next_m_node->chr_index + this->l_fragments)
						{
							search_done = true;
							n_cur_frag_copies++;
						}
						else if(cur_next_r_node->chr_index > cur_next_m_node->chr_index + this->l_fragments)
						{
							// Right node just got out of range of search, it is not necessary to search any more.
							search_done = true;						
						}

						cur_next_r_node = cur_next_r_node->next;
					}
				}
				// If middle node got out of range of possible matching, do not search any more.
				else if(cur_next_m_node->chr_index > cur_l_entry->chr_index + this->l_fragments)
				{
					// Right node just got out of range of search, it is not necessary to search any more.
					search_done = true;
				}
				cur_next_m_node = cur_next_m_node->next;
			}

			// Move to the next left entry.
			cur_l_entry = cur_l_entry->next;
		}
	
		//printf("[%d]         \n", n_cur_frag_copies);
	
		// Dump the count for this position in the mappability file. Note that most probably this will not be needed
		// later, only thing that is necessary is whether this fragment is uniquely mappable or not.
		fwrite(&n_cur_frag_copies, sizeof(int), 1, f_mappability_map);

		base_i++; // Update the position.
    }  //base_i loop.

	fclose(f_mappability_map);
}



