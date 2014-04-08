#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "chromosome.h"

bool sget_line(char* string, int& i, char* buff);

t_chromosome::t_chromosome(char* chromosome_fp)
{
	this->set_quaternary_codes();

	// Add the chromosome in the arguments.
	this->load_chromosome(chromosome_fp);
}

void t_chromosome::set_quaternary_codes()
{
	this->quaternary_codes_per_ASCII_char = new char[256];

	// Init all the values to value of NUC_X, unknown nucleotide.
	memset(this->quaternary_codes_per_ASCII_char, NUC_X, 256 * sizeof(char));

	// Set the values for the 4 nucleotides.
	this->quaternary_codes_per_ASCII_char[(int)'A'] = NUC_A;
    this->quaternary_codes_per_ASCII_char[(int)'a'] = NUC_A;
    this->quaternary_codes_per_ASCII_char[(int)'C'] = NUC_C;
    this->quaternary_codes_per_ASCII_char[(int)'c'] = NUC_C;
    this->quaternary_codes_per_ASCII_char[(int)'G'] = NUC_G;
    this->quaternary_codes_per_ASCII_char[(int)'g'] = NUC_G;
    this->quaternary_codes_per_ASCII_char[(int)'T'] = NUC_T;
    this->quaternary_codes_per_ASCII_char[(int)'t'] = NUC_T;
}

void t_chromosome::load_chromosome(char* chromosome_fp)
{
	printf("Adding chromosome %s\n", chromosome_fp);

	// Initialize the info for new chromosome.
	this->length = 0;
	this->nucs = NULL;
	this->fp = new char[strlen(chromosome_fp) + 2];
	strcpy(this->fp, chromosome_fp);

	// Open the file, count the nucleotides (for determining the size of nucleotide data buffer), allocate the buffer and fill the buffer 
	// with all the nucleotides in (nucs) member.
	FILE* f_chrom = fopen(chromosome_fp, "rb");

	if(f_chrom == NULL)
	{
		printf("Could not open chromosome file @ %s for reading.\n", chromosome_fp);
		exit(0);
	}

	fseek(f_chrom, SEEK_SET, SEEK_END);
	int file_size = (int)(ftell(f_chrom)); // Does not work for chromosomes longer than 4 gigabases.
            fseek(f_chrom, SEEK_SET, SEEK_SET); // Go back to the beginning.
	
	char* file_contents = new char[file_size+10];
	int i_byte = 0;

	// Following makes sure that whole chromosome file is read to memory.
	if(fread(file_contents, file_size, 1, f_chrom) != 1)
	{
		printf("Could load %s into memory.\n", chromosome_fp);
		exit(0);
	}

	file_contents[file_size] = 0; // Treats file contents as a string.

	// Loaded the file to the buffer, close the file.
	fclose(f_chrom);

	// Process file contents.

	// Generate the binary nucleotide data for this chromosome for random access.
	char line_buffer[1000];
	int n_nucs = 0;
	int i_nuc = 0;
	while(sget_line(file_contents, i_nuc, line_buffer))
	{
		//fprintf(stderr, "cur_line:\n%s\n", line_buffer);
		//getc(stdin);

		// Process current line.
		if(strlen(line_buffer) == 0)
		{
		}
		else if(line_buffer[0] == '>')
		{
			// Comment line.
		}
		else
		{
			// Nucleotide data, read.
			for(int i = 0; i < strlen(line_buffer); i++)
			{
				if(line_buffer[i] == 'A' || line_buffer[i] == 'a' ||
					line_buffer[i] == 'C' || line_buffer[i] == 'c' ||
					line_buffer[i] == 'G' || line_buffer[i] == 'g' ||
					line_buffer[i] == 'T' || line_buffer[i] == 't' ||
					line_buffer[i] == 'N' || line_buffer[i] == 'n')
				{
					n_nucs++;
				}
				else if((line_buffer[i] > 'A' && line_buffer[i] < 'Z') ||
					(line_buffer[i] > 'a' && line_buffer[i] < 'z'))
				{
					printf("Adding '%c' in %s\n", line_buffer[i], chromosome_fp);
					n_nucs++;
				}

			} // i loop for current line.
		}
	} // file reading loop.


	// Finished counting the nucleotides, can allocate the buffer to store the nucleotides.

	// Allocate the nucleotides.
	this->nucs = new char[n_nucs + 5];
	printf("Found %d nucleotides.\n", n_nucs);
	this->length = n_nucs;

	// Generate the binary nucleotide data for this chromosome for random access.
	n_nucs = 0; // The chromosome is 1 based.
	this->nucs[n_nucs] = '#'; // The chromosome is 1 based-indexed.
	n_nucs++;
	i_nuc = 0;
            while(sget_line(file_contents, i_nuc, line_buffer))
            {
                    //fprintf(stderr, "cur_line:\n%s\n", line_buffer);
                    //getc(stdin);

                    // Process current line.
                    if(strlen(line_buffer) == 0)
                    {
                    }
                    else if(line_buffer[0] == '>')
                    {
                            // Comment line.
                    }
                    else
                    {
			// Nucleotide data, read.
			for(int i = 0; i < strlen(line_buffer); i++)
			{
                                    if(line_buffer[i] == 'A' || line_buffer[i] == 'a' ||
                                            line_buffer[i] == 'C' || line_buffer[i] == 'c' ||
                                            line_buffer[i] == 'G' || line_buffer[i] == 'g' ||
                                            line_buffer[i] == 'T' || line_buffer[i] == 't' ||
                                            line_buffer[i] == 'N' || line_buffer[i] == 'n')
                                    {
					this->nucs[n_nucs] = line_buffer[i];
                                            n_nucs++;
                                    }
                                    else if((line_buffer[i] > 'A' && line_buffer[i] < 'Z') ||
                                            (line_buffer[i] > 'a' && line_buffer[i] < 'z'))
                                    {
                                            printf("Adding '%c' in %s\n", line_buffer[i], chromosome_fp);
                                            this->nucs[n_nucs] = line_buffer[i];
                                            n_nucs++;
                                    }
			} // i loop for current line.
		}
	} // file reading loop.

	// Dump first 20000 nucleotides, for checking.
/*
	for(int i = 0; i < 20000; i++)
	{
		printf("%c", this->nucs[i]);
	} // i loop
*/
	// Free memory for file contents.
	delete [] file_contents;
}

t_chromosome::~t_chromosome()
{
	// Close the file if working with disk_mapped flag on.
	delete [] this->nucs;
}

// Read a line in the string: Note that this function updates the index on the string.
bool sget_line(char* string, int& i, char* buff)
{
	int i_buff = 0;
	bool read = false;
	while(string[i] != 0 && string[i] != '\n')
	{
		read = true;
		buff[i_buff] = string[i];
		i_buff++;
		i++;
	}

	// Skip new line character.
	if(string[i] == '\n')
	{
		i++;
	}

	// End buffer string.
	buff[i_buff] = 0;

	// If this is the end of the string, return false, since nothing is read.
	if(string[i] == 0 && !read)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

char t_chromosome::x_nuc(int i_nuc)
{
	return(this->nucs[i_nuc]);
}

/*
i_nuc is the 1 based index with respect to the 5' side of reverse strand.
*/
char t_chromosome::x_nuc_rev(int i_nuc)
{
	// Addition of 1 is necessary to compensate for 1 based indexing.
	return(this->complement(this->nucs[this->length + 1 - i_nuc]));
}

char t_chromosome::complement(char nuc)
{
	// Be careful about the following string, if the 
	char complement_nucs[] = "NACGTN";

	// Uses a simple trick to determine the complement: The quaternary coding is symmetric with respect to the complementarity.
	//printf("Nuc is %c (%d)\n", nuc, this->quaternary_codes[(int)nuc]);
	return(complement_nucs[5 - this->quaternary_codes_per_ASCII_char[(int)nuc]]);
}

int t_chromosome::l_chromosome()
{
	return(this->length);
}

// Following functions are for saving disk space when working with disk_mapped flag on.
char t_chromosome::pack_nucs_to_bin(char* nucs, int n_nucs)
{
	char packed_nucs = 0;
	int n_bits_per_nuc = (8 / n_nucs);

	if((n_bits_per_nuc * n_nucs) != 8)
	{
		printf("Warning: The byte division is not exact while packing the nucleotides, n_bits_per_nuc = %d @ %s(%d)\n", n_bits_per_nuc);
		getc(stdin);
	}

	for(int i_nuc = 0; i_nuc < n_nucs; i_nuc++)
	{
		packed_nucs += nucs[i_nuc];

		// If this was not the last nucleotide, shift the byte to left by n_bits_per_nuc bits.
		if(i_nuc != n_nucs-1)
		{
			packed_nucs = packed_nucs << n_bits_per_nuc;
		}
	} // i_nuc loop.

	return(packed_nucs);
}

// Unpack n_nucs many nucleotides from byte.
void t_chromosome::unpack_nucs_from_byte(char byte, char* nucs, int n_nucs)
{
	int n_bits_per_nuc = (8 / n_nucs);

	if((n_bits_per_nuc * n_nucs) != 8)
	{
		printf("Warning: The byte division is not exact while unpacking the nucleotides, n_bits_per_nuc = %d @ %s(%d)\n", n_bits_per_nuc);
		getc(stdin);
	}

	// Generate the mask to extract (n_bits_per_nuc) many least significant bits.
	char nuc_mask = 0;

	// Set up nucleotide mask: Add (n_bits_per_nuc) many 1s at the end of nuc_mask.
	for(int i_bit = 0; i_bit < n_bits_per_nuc; i_bit++)
	{
		nuc_mask += (1 << i_bit);
	} // i_bit loop.

	char packed_nucs = byte;

	for(int i_nuc = 0; i_nuc < n_nucs; i_nuc++)
	{
		nucs[i_nuc] = packed_nucs & nuc_mask;

		// Get rid of last (n_bits_per_nuc) many bits.
		packed_nucs = packed_nucs >> n_bits_per_nuc;
	} // i_nuc loop.
}

bool t_chromosome::fill_first_valid_quaternary_coded_fragment_on_fore_strand(char* fragment, int l_fragment, int& base_i)
{
	while(l_fragment <= (this->length - base_i + 1) &&
		!this->fill_quaternary_coded_fragment_on_fore_strand(fragment, l_fragment, base_i))
	{
		base_i++;
	}

	if(l_fragment <= (this->length - base_i + 1))
	{
		return(true);
	}
	else
	{
		// The filled buffer is not valid, do not use it.
		return(false);
	}
}

// Fill a fragment buffer (of coded nucleotides), using the chromosome information.
bool t_chromosome::fill_quaternary_coded_fragment_on_fore_strand(char* fragment, int l_fragment, int base_i)
{
	// Fill up the fragment.
	for(int i = 0; i < l_fragment; i++)
	{
		fragment[i] = this->quaternary_codes_per_ASCII_char[this->x_nuc(base_i + i)];

		if(fragment[i] == NUC_X || fragment[i] == NUC_COMPX)
		{
			// Could not fill.
			return(false);
		}
	} // i loop

	return(true);
}

/*
base_i is 1 based index with respect to the 5' side of the forward strand.
*/
bool t_chromosome::fill_first_valid_quaternary_coded_fragment_on_rev_strand(char* quad_coded_fragment, int l_fragment, int& base_i)
{
	//while(l_fragment <= (this->length - base_i + 1) &&

	// base_i is decremented at each check.
	while(base_i >= 1 &&
			!this->fill_quaternary_coded_fragment_on_rev_strand(quad_coded_fragment, l_fragment, base_i))
	{
		base_i--;
	}

	if(base_i >= 1)
	{
		return(true);
	}
	else
	{
		// The filled buffer is not valid, do not use it.
		return(false);
	}
}

/* 
Fill a fragment buffer (of coded nucleotides), using the chromosome information. 
base_i is the index on the forward strand that points to the end of the desired fragment.

		base+i
          |
>--------->>>>>>>>>>>>>>>>>>>>>>>---->
<---------<<<<<<<<<<<<<<<<<<<<<<<----<
		  3'<-Desired fragment->5'

The returned fragment is the reverse complement of the forward strand, which is the easiest way
to get the desired fragment.
*/
bool t_chromosome::fill_quaternary_coded_fragment_on_rev_strand(char* quad_coded_fragment, int l_fragment, int base_i)
{
	// Fill up the fragment on the forward strand.
	for(int i = 0; i < l_fragment; i++)
	{
		quad_coded_fragment[i] = this->quaternary_codes_per_ASCII_char[this->x_nuc(base_i + i)];

		// If this is an unknown nucleotide, return 
		if(quad_coded_fragment[i] == NUC_X || quad_coded_fragment[i] == NUC_COMPX)
		{
			// Could not fill.
			return(false);
		}
	}

	// Reverse complement quad_coded_fragment.
	for(int i = 0; i < l_fragment / 2; i++)
	{
		char temp_nuc_code = quad_coded_fragment[i];
		quad_coded_fragment[i] = 5 - quad_coded_fragment[(l_fragment - 1) - i]; // Complement the nucleotide.
		quad_coded_fragment[(l_fragment - 1) - i] = 5 - temp_nuc_code; // Complement the nucleotide.
	} // reverse complementing loop.

	return(true);
}



