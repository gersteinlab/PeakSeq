#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "annot_region_tools.h"
#include "rng.h"
#include "seed_manager.h"
#include "utils.h"
#include "fragment_simulation.h"

const double PI = 3.14159265;

void get_gaussian_cdf_per_region(int l_region, double sigma, vector<double>* cdf)
{
	// Compute the gaussian and sum up to get the cdf.
	//vector<double>* cdf_vals = new vector<double>();
	cdf->push_back(0.0);
	for(int i = 0; i < l_region; i++)
	{
		double coord = (double)i - (((double)l_region)/2);
		double denom = 1 / (pow(2 * PI, .5) * sigma);
		double cur_prob = denom * exp(-1 * coord * coord / (2 * sigma * sigma));

		double cur_cumul_prob = cdf->back() + cur_prob;

		cdf->push_back(cur_cumul_prob);
	} // i loop.

	fprintf(stderr, "CDF end for gaussian distribution is %lf\n", cdf->back());
}

void get_uniform_cdf_per_region(int l_region, vector<double>* cdf)
{
	// Compute the gaussian and sum up to get the cdf.
	//vector<double>* cdf_vals = new vector<double>();
	cdf->push_back(0.0);
	for(int i = 0; i < l_region; i++)
	{
		double cur_prob = (double)i / l_region;

		double cur_cumul_prob = cdf->back() + cur_prob;

		cdf->push_back(cur_cumul_prob);
	} // i loop.

	fprintf(stderr, "CDF end for uniform distribution is %lf\n", cdf->back());
}

void generate_simulated_fragments(char* bed_fp, int n_chip_frags_per_region, int n_ip_frags_per_region, int l_frag, char* chip_op_tagalign_fp, char* ip_op_tagalign_fp)
{
	// Load the regions.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	vector<t_annot_region*>* regions = load_BED(bed_fp);

	FILE* f_chip_tagalign = open_f(chip_op_tagalign_fp, "w");
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		int l_cur_reg = regions->at(i_reg)->end - regions->at(i_reg)->start;
		vector<double>* cdf = new vector<double>();
		get_gaussian_cdf_per_region(l_cur_reg, l_cur_reg / 40, cdf);
		for(int i_frag = 0; i_frag < n_chip_frags_per_region; i_frag++)
		{
			// Generate random indices.
			double cur_rand = rng->random_double_ran3();

			// Find the first index for which the probability is above cdf.
			int i_ind = 0;
			for(int i_cdf = 0; i_cdf < (int)cdf->size(); i_cdf++)
			{
				if(cdf->at(i_cdf) > cur_rand)
				{
					i_ind = i_cdf;
					break;
				}
			} // i_cdf loop.
			
			int genome_coord = regions->at(i_reg)->start + i_ind;

			// Choose the strand randomly.
			char strand = (rng->random_double_ran3() > 0.5)?('+'):('-');

			// Dump the read.
			fprintf(f_chip_tagalign, "%s\t%d\t%d\t.\t1000\t%c\n", regions->at(i_reg)->chrom, genome_coord, genome_coord+l_frag, strand);
		} // i_frag loop.
	} // i_reg loop.
	fclose(f_chip_tagalign);

	// For each position of binding, generate some reads in the input, which is uniformly distributed over the region.
	FILE* f_ip_tagalign = open_f(ip_op_tagalign_fp, "w");
	for(int i_reg = 0; i_reg < (int)regions->size(); i_reg++)
	{
		int l_cur_reg = regions->at(i_reg)->end - regions->at(i_reg)->start;
		vector<double>* cdf = new vector<double>();
		get_uniform_cdf_per_region(l_cur_reg, cdf);
		for(int i_frag = 0; i_frag < n_ip_frags_per_region; i_frag++)
		{
			// Generate random indices.
			double cur_rand = rng->random_double_ran3();

			// Find the first index for which the probability is above cdf.
			int i_ind = 0;
			for(int i_cdf = 0; i_cdf < (int)cdf->size(); i_cdf++)
			{
				if(cdf->at(i_cdf) > cur_rand)
				{
					i_ind = i_cdf;
					break;
				}
			} // i_cdf loop.
			
			int genome_coord = regions->at(i_reg)->start + i_ind;

			// Choose the strand randomly.
			char strand = (rng->random_double_ran3() > 0.5)?('+'):('-');

			// Dump the read.
			fprintf(f_ip_tagalign, "%s\t%d\t%d\t.\t1000\t%c\n", regions->at(i_reg)->chrom, genome_coord, genome_coord+l_frag, strand);
		} // i_frag loop.
	} // i_reg loop.
	fclose(f_ip_tagalign);
}


