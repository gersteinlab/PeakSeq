#include <stdio.h>
#include <stdlib.h>
#include "mapped_read_tools.h"
#include "input_normalization.h"
#include "enrichment_profile.h"
#include "simulated_background.h"
#include "config.h"
#include "compare_signal_tracks.h"
#include "peakseq.h"
#include "fragment_simulation.h"
#include <time.h>
#include "utils.h"
#include <string.h>

#define MAX_BUFF_SIZE (10000)

int main(int argc, char* argv[])
{
	// Load chromosome file list.
        if(argc < 3)
        {
                printf("\n\
USAGE: %s [option] [data files]\n \
	Options:\n\
	-preprocess [Input format ('SAM'/'ELAND'/'tagAlign'/'bowtie')] [Path to the file that contains mapped tag data ('stdin' for piping)] [output directory for parsed reads]\n\
	-peak_select [Configuration file path]\n\
Configuration file example: \n\
	chromosome_list_file [path to the chromosome list file]\n\
	mappability_map_file [path to the mappability map file]\n\
	minimum_gap_between_peaks [minimum gap between the peaks] \n\
	extended_region_size [extended region size for correction of # of tags for input experiments] \n\n\
** Important notes about chromosome ids for -peak_select option:\n\
	1. The ids in the mappability map file must be consistent with the ids in the id list file.\n\
	2. The ids in the id list file must correspond to the suffixes of the binary mapped tag data files. This is how PeakSeq reads the correct tag data for a chromosome.\n", argv[0]);

                exit(0);
        }

	clock_t start = clock();

	// Other preprocessing/input validation options go here.
	if(strcmp(argv[1], "-peak_select") == 0)
	{
		if(argc != 3)
		{
			printf("USAGE: -peak_select [Configuration file path]\n");
			exit(0);
		}

		t_config* config = new t_config(argv[2]);
		char chr_list_fp[__MAX_PATH];

		int enrichment_mapped_fragment_length;
        if(!config->get_int_val("Enrichment_mapped_fragment_length", enrichment_mapped_fragment_length))
        {
                printf("Could not find Enrichment_mapped_fragment_length entry in the config file %s\n", argv[1]);
                exit(0);
        }

		double target_FDR;
        if(!config->get_double_val("target_FDR", target_FDR))
        {
                printf("Could not find target_FDR entry in the config file %s\n", argv[1]);
                exit(0);
        }

		double max_Qvalue;
        if(!config->get_double_val("max_Qvalue", max_Qvalue))
        {
                printf("Could not find max_Qvalue entry in the config file %s\n", argv[1]);
                exit(0);
        }

		int min_thresh = 1;
		int max_thresh = 1000;
		int n_sims = 10;
        if(!config->get_int_val("N_Simulations", n_sims))
        {
                printf("Could not find N_Simulations entry in the config file %s\n", argv[1]);
                exit(0);
        }

		int min_gap_per_bw_peaks;
        if(!config->get_int_val("Minimum_interpeak_distance", min_gap_per_bw_peaks))
        {
                printf("Could not find Minimum_interpeak_distance entry in the config file %s\n", argv[1]);
                exit(0);
        }

		double P_f = 0.0f; // Percentage to exclude from the scaling factor estimation for input data.

        char mappability_map_fp[__MAX_PATH];
        if(!config->get_str_val("Mappability_map_file", mappability_map_fp))
        {
                printf("Could not find Mappability_map_file entry in the config file %s\n", argv[1]);
                exit(0);
        }

		// Read the background model.
		char background_model[100];
		if(!config->get_str_val("Background_model", background_model))
		{
            printf("Could not find Background_model entry in the config file %s\n", argv[1]);
            exit(0);				
		}
		else
		{
			if(strcmp(background_model, "Poisson") != 0 &&
				strcmp(background_model, "Simulated") != 0)
			{
				printf("Use 'Poisson' or 'Simulated' for Background_model option.\n");
				exit(0);
			}
		}

		vector<char*>* chip_seq_reads_data_dirs = config->get_val_list("ChIP_Seq_reads_data_dirs");
        if(chip_seq_reads_data_dirs == NULL)
        {
                printf("Could not find ChIP_Seq_reads_data_dirs entry in the config file %s\n", argv[1]);
                exit(0);
        }		

		// Get the chromosoe id's.
		sprintf(chr_list_fp, "%s/chr_ids.txt", chip_seq_reads_data_dirs->at(0));

		vector<char*>* input_reads_data_dirs = config->get_val_list("Input_reads_data_dirs");
        if(input_reads_data_dirs == NULL)
        {
                printf("Could not find Input_reads_data_dirs entry in the config file %s\n", argv[1]);
                exit(0);
        }

		char exp_id[MAX_BUFF_SIZE];
		if(!config->get_str_val("Experiment_id", exp_id))
		{
			printf("Could not find Experiment_id entry in the config file %s\n", argv[1]);
			exit(0);
		}

		char peak_profs_op_dir[__MAX_PATH];
		peak_profs_op_dir[0] = 0;
                if(!config->get_str_val("Peak_profile_output_dir", peak_profs_op_dir))
                {
                        //printf("Could not find Peak_profile_output_dir in the config file %s\n", argv[1]);
                }

		char narrowPeak_op_fp[__MAX_PATH];
		narrowPeak_op_fp[0] = 0;
                if(!config->get_str_val("narrowPeak_output_file_path", narrowPeak_op_fp))
                {
                        //printf("Could not find narrowPeak_output_file_path in the config file %s\n", argv[1]);
                }

               	int rand_seed = -1; 
                if(!config->get_int_val("Simulation_seed", rand_seed))
                {
                        //printf("Could not find Simulation_seed in the config file %s\n", argv[1]);
                        //exit(0);
                }
				// Generates the read file paths from chr id's.
				peakseq(exp_id,
						chr_list_fp,
		        			enrichment_mapped_fragment_length,
							background_model,
							target_FDR,
							max_Qvalue,
							min_thresh,
							max_thresh,
							n_sims,
							min_gap_per_bw_peaks,
							P_f,
							chip_seq_reads_data_dirs,
							input_reads_data_dirs,
							mappability_map_fp,
							peak_profs_op_dir,
							narrowPeak_op_fp,
							rand_seed); // Checks if the peak profiles op dir is 
	}
	else if(strcmp(argv[1], "-preprocess") == 0)
	{
		if(argc != 5)
		{
			printf("USAGE: %s -preprocess [Input format ('SAM'/'ELAND'/'tagAlign'/'bowtie')] [Path to the file that contains mapped tag data ('stdin' for piping)] [output directory for parsed reads]\n", argv[0]);
			exit(0);
		}
/*-preprocess 
	[Input format ('SAM'/'ELAND')] 
	[Path for the file that contains chromosome id list] 
	[Path to the file that contains mapped tag data ('stdin' for piping)] 
	[output directory for parsed reads]
*/
		char ip_format[MAX_BUFF_SIZE];
		strcpy(ip_format, argv[2]);

		char chip_seq_eland_op_fp[__MAX_PATH];
		strcpy(chip_seq_eland_op_fp, argv[3]);
	
		char parsed_reads_op_dir[__MAX_PATH];
		strcpy(parsed_reads_op_dir, argv[4]);

		printf("Preprocessing:\ninput format: %s\nchip_seq_eland_op_fp: %s\nparsed_reads_op_dir: %s\n", 
			ip_format, chip_seq_eland_op_fp, parsed_reads_op_dir);

		// Do not do validation while dumping the binary read files.
		if(strcmp(ip_format, "ELAND") == 0)
		{
			// Read the ELAND file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_ELAND_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(chip_seq_eland_op_fp, parsed_reads_op_dir, preprocess_ELAND_read_line, false);
		}
		else if(strcmp(ip_format, "SAM") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_SAM_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(chip_seq_eland_op_fp, parsed_reads_op_dir, preprocess_SAM_read_line, false);
		}
		else if(strcmp(ip_format, "tagAlign") == 0)
		{
			// Read the SAM file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_tagAlign_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(chip_seq_eland_op_fp, parsed_reads_op_dir, preprocess_tagAlign_read_line, false);
		}
		else if(strcmp(ip_format, "bowtie") == 0)
		{
			// Read the bowtie file specified at the command line, separate it: Generate a mapped _reads file for each chromosome.
			//parse_bowtie_formatted_mapped_reads_file(chr_fps, parsed_reads_op_dir, chip_seq_eland_op_fp);
			preprocess_mapped_reads_file(chip_seq_eland_op_fp, parsed_reads_op_dir, preprocess_bowtie_read_line, false);
		}
		else
		{
			printf("Unknown format for the mapped read file name, use SAM/ELAND.\n");
			exit(0);
		}
	}
	else if(strcmp(argv[1], "-simulate_fragments") == 0)
	{
		if(argc != 8)
		{
			printf("USAGE: -simulate_fragments [Path for bed file where the reads will be localized to] [# of fragments to per region] [fragment length] [ChIP output file path] [Input output file path]\n");
			exit(0);
		}

		char* bed_fp = argv[2];
		int n_chip_frags_per_region = atoi(argv[3]);
		int n_ip_frags_per_region = atoi(argv[4]);
		int l_frag = atoi(argv[5]);
		char* chip_op_tagalign_fp = argv[6];
		char* ip_op_tagalign_fp = argv[7];
		generate_simulated_fragments(bed_fp, n_chip_frags_per_region, n_ip_frags_per_region, l_frag, chip_op_tagalign_fp, ip_op_tagalign_fp);
	}
	//else if(strcmp(argv[1], "-validate_fragments") == 0)
	//{
 //       char chr_list_fp[__MAX_PATH];
 //       strcpy(chr_list_fp, argv[2]);

	//	char parsed_reads_dir[__MAX_PATH];
	//	strcpy(parsed_reads_dir, argv[3]);

 //       // Load the list of mappability map file paths.
 //       vector<char*>* chr_fps = new vector<char*>();
 //       FILE* f_chr_list = open_f(chr_list_fp, "r");
 //       if(f_chr_list == NULL)
 //       {
 //               printf("Could not open chrososome list file %s\n", chr_list_fp);
 //               exit(0);
 //       }

 //       char cur_fp[__MAX_PATH];
 //       while(fscanf(f_chr_list, "%s", cur_fp) == 1)
 //       {
 //               char* new_chr_fp = new char[strlen(cur_fp) + 4];
 //               sprintf(new_chr_fp, cur_fp);
 //               chr_fps->push_back(new_chr_fp);
 //               //printf("%s\n", new_chr_fp);
 //       }
 //       fclose(f_chr_list);

	//	// Request validation.
 //       //validate_dump_binary_files(chr_fps, parsed_reads_dir, true);		
	//}
	else
	{
		printf("Could not understand the option %s, use -peak_select or -preprocess.\n", argv[1]);
		exit(0);
	}

	clock_t end = clock();

	// Dump the beacon file, which indicates successful exit from PeakSeq.
	FILE* f_beacon = open_f("beacon.log", "w");
	fprintf(f_beacon, "Finished in %lf seconds.\n", ((double)end - (double)start) / CLOCKS_PER_SEC);
	fclose(f_beacon);

	return(0);
}


