#include <stdio.h>
#include <stdlib.h>
#include "peakseq_output.h"
#include <string.h>
#include "utils.h"

// Define the static object.
t_peakseq_output* t_peakseq_output::singleton = NULL;

/*
class t_peak_seq_output
{
        // Check if the object is allocated, if so do not do anything. Otherwise setup the data: Initialize the output file, copy the id.
        // Set the current chromosome id to emptry string.
        t_peakseq_output(char* _chip_seq_experiment_id);
        ~t_peakseq_output();

        static t_output* singleton;
        char cur_chr_id[100];
        char chip_seq_experiment_id[100];

        static void set_cur_chr_id(char* new_chr_id);
        static void dump(char* result_string); // Dump a result string.
}
*/

t_peakseq_output::t_peakseq_output(int i)
{
	this->cur_chr_id = NULL;
	this->chip_seq_experiment_id = NULL;
	//memset(this->cur_chr_id, 0, 100);
	//memset(this->chip_seq_experiment_id, 0, 100);
}

t_peakseq_output::t_peakseq_output(char* _chip_seq_experiment_id)
{
	if(t_peakseq_output::singleton == NULL)
	{
		//printf("Initializing peakseq output object with %s\n", _chip_seq_experiment_id);
		t_peakseq_output::singleton = new t_peakseq_output(0); // Cannot use the same constructor.
		t_peakseq_output::singleton->chip_seq_experiment_id = new char[strlen(_chip_seq_experiment_id) + 3];
		strcpy(t_peakseq_output::singleton->chip_seq_experiment_id, _chip_seq_experiment_id);

        //char* op_fp = new char[strlen(t_peakseq_output::singleton->chip_seq_experiment_id) + 10];
		char op_fp[__MAX_PATH];
        sprintf(op_fp, "%s.txt", t_peakseq_output::singleton->chip_seq_experiment_id);		
        FILE* f_op = open_f(op_fp, "w");
		fclose(f_op);

		t_peakseq_output::singleton->narrowPeak_op_fp = new char[strlen(t_peakseq_output::singleton->chip_seq_experiment_id) + strlen("_narrowPeak.txt") + 3];
		sprintf(t_peakseq_output::singleton->narrowPeak_op_fp, "%s_narrowPeak.txt", t_peakseq_output::singleton->chip_seq_experiment_id);
		FILE* f_narrowPeak_op = open_f(t_peakseq_output::singleton->narrowPeak_op_fp, "w");
		fclose(f_narrowPeak_op);
		//printf("Initialized peakseq output object.\n");
	}
}

/*
Constructor with narrowPeak output file path.
*/
t_peakseq_output::t_peakseq_output(char* _chip_seq_experiment_id, char* _narrowPeak_op_fp)
{
	if(t_peakseq_output::singleton == NULL)
	{
		//printf("Initializing peakseq output object with %s\n", _chip_seq_experiment_id);
		t_peakseq_output::singleton = new t_peakseq_output(0); // Cannot use the same constructor.
		t_peakseq_output::singleton->chip_seq_experiment_id = new char[strlen(_chip_seq_experiment_id) + 3];
		strcpy(t_peakseq_output::singleton->chip_seq_experiment_id, _chip_seq_experiment_id);

	        char op_fp[__MAX_PATH];
	        sprintf(op_fp, "%s.txt", t_peakseq_output::singleton->chip_seq_experiment_id);		
	        FILE* f_op = open_f(op_fp, "w");
		fclose(f_op);

		//char narrow_peak_op_fp[500];
		//sprintf(narrow_peak_op_fp, "%s_narrowPeak.txt", t_peakseq_output::singleton->chip_seq_experiment_id);
		t_peakseq_output::singleton->narrowPeak_op_fp = new char[strlen(_narrowPeak_op_fp) + 3];
                strcpy(t_peakseq_output::singleton->narrowPeak_op_fp, _narrowPeak_op_fp);
                FILE* f_narrowPeak_op = open_f(t_peakseq_output::singleton->narrowPeak_op_fp, "w");
                fclose(f_narrowPeak_op);
	}
}

void t_peakseq_output::set_cur_chr_id(char* new_chr_id)
{
	if(t_peakseq_output::singleton == NULL)
	{
		printf("The output object is not initialized.\n");
		exit(0);
	}
	else
	{
		t_peakseq_output::singleton->cur_chr_id = new char[strlen(new_chr_id) + 3];
		strcpy(t_peakseq_output::singleton->cur_chr_id, new_chr_id);
	}
}

void t_peakseq_output::dump_PeakSeq_title(char* title_string)
{
        if(t_peakseq_output::singleton == NULL)
        {
                printf("The output object is not initialized.\n");
                exit(0);
        }
        else
        {
                //strcpy(t_peakseq_output::singleton->cur_chr_id, new_chr_id);
                char op_fp[__MAX_PATH];
                sprintf(op_fp, "%s.txt", t_peakseq_output::singleton->chip_seq_experiment_id);
                FILE* f_op = open_f(op_fp, "a");
                fprintf(f_op, "%-20s\n", title_string);
                fclose(f_op);
        }
}

void t_peakseq_output::dump_PeakSeq_line(char* result_string)
{
    if(t_peakseq_output::singleton == NULL)
    {
            printf("The output object is not initialized.\n");
            exit(0);
    }
    else
    {
		//strcpy(t_peakseq_output::singleton->cur_chr_id, new_chr_id);
		char op_fp[__MAX_PATH];
		sprintf(op_fp, "%s.txt", t_peakseq_output::singleton->chip_seq_experiment_id);
		FILE* f_op = open_f(op_fp, "a");
		fprintf(f_op, "%s\n", result_string);
		fclose(f_op);
    }
}

void t_peakseq_output::dump_narrowPeak_line(char* result_string)
{
    if(t_peakseq_output::singleton == NULL)
    {
            printf("The output object is not initialized.\n");
            exit(0);
    }
    else
    {
		//sprintf(narrow_peak_op_fp, "%s_narrowPeak.txt", t_peakseq_output::singleton->chip_seq_experiment_id);
		FILE* f_narrowPeak_op = open_f(t_peakseq_output::singleton->narrowPeak_op_fp, "a");
		//fprintf(f_narrowPeak_op, "%-15s %s\n", t_peakseq_output::singleton->cur_chr_id, result_string);
		fprintf(f_narrowPeak_op, "%s\n", result_string);
		fclose(f_narrowPeak_op);
    }
}


t_peakseq_output::~t_peakseq_output()
{
}


