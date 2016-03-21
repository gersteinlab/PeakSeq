#include <stdio.h>
#include <stdlib.h>
#include "enrichment_profile.h"
#include "mapped_read_tools.h"
#include <math.h>
#include "peak_region.h"
#include <algorithm>
#include "chip_seq_chr_data.h"
#include "utils.h"
#include "ansi_string.h"
#include <string.h>

void t_enrichment_profile::add_fragment_ends(int start, int length)
{
	if(this->frag_ends == NULL)
	{
		this->frag_ends = new vector<t_frag_end*>();
	}

	t_frag_end* frag_start = (t_frag_end*)malloc(sizeof(t_frag_end));
	frag_start->pos = start;
	frag_start->side = FRAG_START;

	t_frag_end* frag_end = (t_frag_end*)malloc(sizeof(t_frag_end));
	frag_end->pos = start + length; // Note that this position points to the nucleotide where profile height goes down. Therefore this is 1 nucleotide after the end of the fragment. The fragment ends at (start + length - 1).
	frag_end->side = FRAG_END;

	// Add these ends.
        this->frag_ends->push_back(frag_start);
	this->frag_ends->push_back(frag_end);
}

//void t_enrichment_profile::get_min_max_height(int& min_height, int& max_height, int min_i_nuc, int max_i_nuc)
//{
//	min_height = 2000;
//	max_height = 0;
//	int height = 0;
//        for(int i_frag = 0; i_frag < this->frag_ends->size(); i_frag++)
//        {
//		height += this->frag_ends->at(i_frag)->side;
//
//		if(this->frag_ends->at(i_frag)->pos >= min_i_nuc && 
//			this->frag_ends->at(i_frag)->pos < max_i_nuc)
//		{
//			if(height > max_height)
//			{
//				max_height = height;
//			}
//	
//			if(height > 0 &&
//				height < min_height)
//			{
//				min_height = height;
//			}
//		}
//	} // i_frag	
//}

/////*
////Dump the profile that spans the requested region.
////*/
//vector<t_profile_site*>* t_enrichment_profile::buffer_profile(int min_i_nuc, int max_i_nuc)
//{
//	int height = 0;
//
//	vector<t_profile_site*>* profile_sites = new vector<t_profile_site*>();
//
//	int prev_pos = 0;
//	int prev_height = 0;
//	int cur_pos = 1;
//    for(int i_end = 0; i_end < this->frag_ends->size(); i_end++)
//    {
//		/*
//		A new region is added between prev_pos and this->frag_ends->at(i_end). If the requested region overlaps
//		with the added region, dump the previous position. Note that the height in the added region is the hiehgt
//		that is updated at the previous position.
//		*/
//        // Update the current height: Note that there can be multiple ends at this position.
//        height += this->frag_ends->at(i_end)->side;
//        cur_pos = this->frag_ends->at(i_end)->pos;
//
//        // Now add all the other ends at this position.
//        while((i_end + 1) < this->frag_ends->size() &&
//                    cur_pos == this->frag_ends->at(i_end+1)->pos)
//        {
//                i_end++;
//                height += this->frag_ends->at(i_end)->side;
//        }
//
//		if(prev_pos >= max_i_nuc)
//		{
//			// The profile passed the requested region completely, can return from the function, dump the previous position and height and return, this is necessary to account for the last position.
//                        //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);
//
//			t_profile_site* new_site = new t_profile_site();
//			new_site->i_nuc = prev_pos;
//			new_site->height = prev_height;
//			profile_sites->push_back(new_site);
//
//			break;
//		}
//		else if(cur_pos < min_i_nuc)
//		{
//			// Still did not reach the requested region.
//		}
//		else
//		{
//			// The newly added region overlaps with the requested region, dump the previous region and height, since 
//			// the region corresponds to the height of the revious pos.
//                        //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);
//			t_profile_site* new_site = new t_profile_site();
//			new_site->i_nuc = prev_pos;
//			new_site->height = prev_height;
//			profile_sites->push_back(new_site);
//		}
//	
//		// Update the previous position before moving to the next region.
//		prev_pos = cur_pos;
//		prev_height = height;
///*
//                if(this->frag_ends->at(i_end)->pos >= min_i_nuc &&
//                        this->frag_ends->at(i_end)->pos < max_i_nuc)
//                {
//			fprintf(f_prof, "%d %d\n", this->frag_ends->at(i_end)->pos, height);
//                }
//*/
//        } // i_frag
//
//		return(profile_sites);
//}

/*
Dump the profile that spans the requested region.
*/
void t_enrichment_profile::dump_profile(char* fp, 
										int min_i_nuc, 
										int max_i_nuc)
{
	int height = 0;
	FILE* f_prof = open_f(fp, "w");
	if(f_prof == NULL)
	{
		printf("Could not open %s for writing profile data.\n", fp);
		exit(0);
	}

	int prev_pos = 0;
	int prev_height = 0;
	int cur_pos = 1;
        for(int i_end = 0; i_end < (int)this->frag_ends->size(); i_end++)
        {
		/*
		A new region is added between prev_pos and this->frag_ends->at(i_end). If the requested region overlaps
		with the added region, dump the previous position. Note that the height in the added region is the hiehgt
		that is updated at the previous position.
		*/
                // Update the current height: Note that there can be multiple ends at this position.
                height += this->frag_ends->at(i_end)->side;
                cur_pos = this->frag_ends->at(i_end)->pos;

                // Now add all the other ends at this position.
                while((i_end + 1) < (int)this->frag_ends->size() &&
                         cur_pos == this->frag_ends->at(i_end+1)->pos)
                {
                        i_end++;
                        height += this->frag_ends->at(i_end)->side;
                }

		if(prev_pos >= max_i_nuc)
		{
			// The profile passed the requested region completely, can return from the function, dump the previous position and height and return, this is necessary to account for the last position.
                        fprintf(f_prof, "%s\t%d\t%d\n", this->chr_id, prev_pos, prev_height);

			break;
		}
		else if(cur_pos < min_i_nuc)
		{
			// Still did not reach the requested region.
		}
		else
		{
			// The newly added region overlaps with the requested region, dump the previous region and height, since 
			// the region corresponds to the height of the revious pos.
                        //fprintf(f_prof, "%d %d\n", prev_pos, prev_height);
						//fprintf(f_prof, "%s\t%d\t%d\n", chr_id, prev_pos, prev_height);
						fprintf(f_prof, "%s\t%d\t%d\n", this->chr_id, prev_pos, prev_height);
		}
	
		// Update the previous position before moving to the next region.
		prev_pos = cur_pos;
		prev_height = height;
/*
                if(this->frag_ends->at(i_end)->pos >= min_i_nuc &&
                        this->frag_ends->at(i_end)->pos < max_i_nuc)
                {
			fprintf(f_prof, "%d %d\n", this->frag_ends->at(i_end)->pos, height);
                }
*/
        } // i_frag
	fclose(f_prof);
}

// This constructor just builds an empty profile to be populated later.
t_enrichment_profile::t_enrichment_profile()
{
	this->frag_ends = NULL;
}

t_enrichment_profile::t_enrichment_profile(char* _chr_id, t_mapped_fragment** fragments, int enrichment_mapped_fragment_length)
{
	this->chr_id = new char[strlen(_chr_id) + 2];
	strcpy(this->chr_id, _chr_id);

        // Initialize counts.
        this->frag_ends = new vector<t_frag_end*>();

        //printf("Building profile with %d fragments.\n", fragments->size());

        // Proces each fragment in the fragments list.
	int i_frag = 0;
        //for(int i_frag = 0; i_frag < fragments->size(); i_frag++)
	while(fragments[i_frag] != NULL)
        {
                t_mapped_fragment* cur_frag = fragments[i_frag];

                if(cur_frag->strand_char == 'F')
                {
						if(enrichment_mapped_fragment_length < cur_frag->sequenced_fragment_length)
						{
							this->add_fragment_ends(cur_frag->base_index, cur_frag->sequenced_fragment_length);
						}
						else
						{
							this->add_fragment_ends(cur_frag->base_index, enrichment_mapped_fragment_length);
						}

                        // Add two points: One for beginning of the fragment and one for the end:
                        //this->add_fragment_ends(cur_frag->base_index - 1, enrichment_mapped_fragment_length);
                        
                }
                else if(cur_frag->strand_char == 'R') // The base indices of reverse strand fragments point to 5' end of the tags.
                {
					if(enrichment_mapped_fragment_length < cur_frag->sequenced_fragment_length)
					{
						this->add_fragment_ends(cur_frag->base_index - (cur_frag->sequenced_fragment_length - 1), cur_frag->sequenced_fragment_length);
					}
					else
					{
						this->add_fragment_ends(cur_frag->base_index - (enrichment_mapped_fragment_length - 1), enrichment_mapped_fragment_length);
					} 
                }
                else
                {
                        printf("Could not resolve strand char '%c' for current fragment, i_frag = %d @ %s(%d)\n", cur_frag->strand_char, i_frag, __FILE__, __LINE__);
                        exit(0);
                }

		i_frag++;
        } // i_frag loop.

        // Sort the fragment ends.
        sort(this->frag_ends->begin(), this->frag_ends->end(), sort_ends);
        //fprintf(stderr, "Done!\n");
}

t_enrichment_profile::t_enrichment_profile(char* _chr_id, vector<t_mapped_fragment*>* fragments, int enrichment_mapped_fragment_length)
{
	this->chr_id = new char[strlen(_chr_id) + 2];
	strcpy(this->chr_id, _chr_id);

        // Initialize counts.
        this->frag_ends = new vector<t_frag_end*>();

	//printf("Building profile with %d fragments.\n", fragments->size());

        // Proces each fragment in the fragments list.
        for(int i_frag = 0; i_frag < (int)fragments->size(); i_frag++)
        {
                t_mapped_fragment* cur_frag = fragments->at(i_frag);

                if(cur_frag->strand_char == 'F')
                {
					// Add two points: One for beginning of the fragment and one for the end:
					//this->add_fragment_ends(cur_frag->base_index - 1, enrichment_mapped_fragment_length);
					if(enrichment_mapped_fragment_length < cur_frag->sequenced_fragment_length)
					{
						this->add_fragment_ends(cur_frag->base_index, cur_frag->sequenced_fragment_length);
					}
					else
					{
						this->add_fragment_ends(cur_frag->base_index, enrichment_mapped_fragment_length);
					}
                }
                else if(cur_frag->strand_char == 'R')
                {
					fprintf(stderr, "Cannot have reverse strand fragment ends, aborting.\n");
					exit(0);
					if(enrichment_mapped_fragment_length < cur_frag->sequenced_fragment_length)
					{
						this->add_fragment_ends(cur_frag->base_index - (cur_frag->sequenced_fragment_length - 1), cur_frag->sequenced_fragment_length);
					}
					else
					{
						this->add_fragment_ends(cur_frag->base_index - (enrichment_mapped_fragment_length - 1), enrichment_mapped_fragment_length);
					} 
                }
                else
                {
                        printf("Could not resolve strand char '%c' for current fragment, i_frag = %d(%ld) @ %s(%d)\n", cur_frag->strand_char, i_frag, fragments->size(), __FILE__, __LINE__);
                        exit(0);
                }
        } // i_frag loop.

        // Sort the fragment ends.
        sort(this->frag_ends->begin(), this->frag_ends->end(), sort_ends);
        //fprintf(stderr, "Done!\n");
}

void delete_frag_ends(vector<t_frag_end*>* frag_ends)
{
	for(int i_end = 0; i_end < (int)frag_ends->size(); i_end++)
	{
		delete(frag_ends->at(i_end));
	} // i_end loop.

	delete(frag_ends);
}

void add_fragment_ends(vector<t_frag_end*>* frag_ends, int start, int length)
{
	t_frag_end* frag_start = new t_frag_end();
	frag_start->pos = start;
	frag_start->side = FRAG_START;

	t_frag_end* frag_end = new t_frag_end();
	frag_end->pos = start + length; // Note that this position points to the nucleotide where profile height goes down. Therefore this is 1 nucleotide after the end of the fragment. The fragment ends at (start + length - 1).
	frag_end->side = FRAG_END;

	// Add these ends.
	frag_ends->push_back(frag_start);
	frag_ends->push_back(frag_end);
}

// Return the fragment ends for a valid set of fragments. Note that the tag length extension is not applied here any more.
vector<t_frag_end*>* get_frag_ends(vector<t_mapped_fragment*>* fragments)
{
        // Initialize counts.
        vector<t_frag_end*>* frag_ends = new vector<t_frag_end*>();

		//printf("Building profile with %d fragments.\n", fragments->size());

        // Proces each fragment in the fragments list.
        for(int i_frag = 0; i_frag < (int)fragments->size(); i_frag++)
        {
				t_mapped_fragment* cur_frag = fragments->at(i_frag);

                if(cur_frag->strand_char == 'F')
                {
					// Add two points: One for beginning of the fragment and one for the end:
					//this->add_fragment_ends(cur_frag->base_index - 1, enrichment_mapped_fragment_length);
					add_fragment_ends(frag_ends, cur_frag->base_index, cur_frag->sequenced_fragment_length);
                }
                else if(cur_frag->strand_char == 'R')
                {
					add_fragment_ends(frag_ends, cur_frag->base_index - (cur_frag->sequenced_fragment_length - 1), cur_frag->sequenced_fragment_length);
                }
                else
                {
                        printf("Could not resolve strand char '%c' for current fragment, i_frag = %d(%ld) @ %s(%d)\n", cur_frag->strand_char, i_frag, fragments->size(), __FILE__, __LINE__);
                        exit(0);
                }
        } // i_frag loop.

        // Sort the fragment ends.
        sort(frag_ends->begin(), frag_ends->end(), sort_ends);
        //fprintf(stderr, "Done!\n");

		return(frag_ends);
}

bool sort_ends(t_frag_end* end1, t_frag_end* end2)
{
	return(end1->pos < end2->pos);
}

void t_enrichment_profile::count_peaks_per_frags(int* n_peaks_per_thresh, vector<int>* frag_begs, vector<int>* frag_ends, int min_thresh, int max_thresh, int min_gap_bw_peaks)
{
        signed int* next_possible_peak_starts = new signed int[max_thresh - min_thresh + 1];
        bool* above_threshold = new bool[max_thresh - min_thresh + 1];
        int* peak_starts = new int[max_thresh - min_thresh + 1];
        for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
        {
                next_possible_peak_starts[thresh - min_thresh] = -10000;
                n_peaks_per_thresh[thresh - min_thresh] = 0;
                above_threshold[thresh - min_thresh] = false;
                peak_starts[thresh - min_thresh ] = -10000;
        } // thresh loop.

    	// Go over the profile and count the peaks.
    	signed int height = 0;
        int cur_pos = 0; // The current fragment end position. This is the position for which the enrichment pfofile is last updated.
	int i_ends = 0; // current index to the fragment begins list.
	int i_begs = 0; // current index to the fragment ends list.
	signed int diff_score = 0;

	while(i_ends < (int)frag_ends->size())
	{
		// Following function handles updating the indices in the list of ends.
		get_height_change(diff_score,
				cur_pos,
				i_ends,
				i_begs,
				frag_begs,
				frag_ends);

		height += diff_score;

/*
		if(height != 0)
		{
			FILE* f_rand_prof = open_f("random_profile.txt", "a");
			fprintf(f_rand_prof, "%d %d\n", cur_pos, height);
			fclose(f_rand_prof);
		}
*/
		// At this point, the height is changed. 
		// For each threhold that is lower than this height, this is either a new peak (if signal was below the threshold) otherwise 
		// must check if the closest peak was at least 200 bps away.
		// If height change is 0, the profile does not change, do not process this case.
		for(int thresh = min_thresh; 
			thresh <= max_thresh; 
			thresh++)
		{
			if(height >= thresh)
			{
				if(above_threshold[thresh - min_thresh])
				{
					// The signal was above the threshold in the last signal update, thus it has to be above this threshold. Do not need to do anything at this point.
				}
				else // Was below the threshold.
				{
					// check if this peak is far enough from the last peak.
					if(cur_pos > next_possible_peak_starts[thresh - min_thresh])
					{
						// This starts a new peak, and ends the previous pesk. Check if the last peak is in the requested window:
						int last_peak_end = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
						int last_peak_start = peak_starts[thresh - min_thresh];
		
						// This check is necessary to skip adding the negative initialized values as a peak.
						if(last_peak_start > 0 && last_peak_end > 0)
						{	
							// If there is overlap between the peak and the investigated region, update the count by one.
							n_peaks_per_thresh[thresh - min_thresh]++;
						}
						//printf("Added %d<%d-%d>\n", thresh, last_peak_start, last_peak_end);

						// Mark the start of a new peak for this threshold.
						peak_starts[thresh - min_thresh] = cur_pos;
					}
				}
	
				// Set the new signal state as above threshold for this peak.
				above_threshold[thresh - min_thresh] = true;
			}
			else
			{			
				// If the signal was above this threshold previously, then this is the end of a peak, mark the next possible peak start for this trheshold.
				if(above_threshold[thresh - min_thresh])
				{
					next_possible_peak_starts[thresh - min_thresh] = cur_pos + min_gap_bw_peaks;
				}

				above_threshold[thresh - min_thresh] = false;
			}
		} // thresh loop
	} // i_end loop.

	// Ensure that the last peak is accounted for:
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		// The signal was below the threshold in the last height update. There are two possibilities: 
		// 1. There was really a peak and algorithm was waiting for another one to check if they need to be emerged.
		// 2. There was never a peak for this threshold.
		if(!above_threshold[thresh - min_thresh])
		{
			int last_peak_end = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
			//int last_peak_start = peak_starts[thresh - min_thresh];

			if(last_peak_end > 0)
				n_peaks_per_thresh[thresh - min_thresh]++;
		}
		else
		{
			// The peak was above threshold in the last position. This means that all the positions remaining after this
			// will be incremented by one.
			n_peaks_per_thresh[thresh - min_thresh]++;
		}

		//printf("%d: %d\n", thresh, n_peaks_per_thresh[thresh - min_thresh]);
	} // i_thresh loop.

//	exit(0);

	//fprintf(stderr, "Ended counting!\n");
	//getc(stdin);

	// Free the array memory.
    delete [] next_possible_peak_starts;
    delete [] above_threshold;
    delete [] peak_starts;
	next_possible_peak_starts = NULL;
	above_threshold = NULL;
	peak_starts = NULL;
}

/*
Following function computes the score for the next index where the height changes. 
*/
void get_height_change(signed int& score, 
			int& index, 
			int& i_ends, 
			int& i_begs, 
			vector<int>* frag_begs, 
			vector<int>* frag_ends)
{
	score = 0;
	
	//printf("Getting height change: begins top: %d, ends top: %d\n", frag_begs->at(i_begs), frag_ends->at(i_ends));
	//getc(stdin);

	if(i_begs < (int)frag_begs->size() &&
		frag_begs->at(i_begs) < frag_ends->at(i_ends))
	{
		index = frag_begs->at(i_begs);
		i_begs++;
		score++;

		// Go over all the ends that have the same index in the begins list.
		while(i_begs < (int)frag_begs->size() &&
			index == frag_begs->at(i_begs))
		{
                        //printf("Begins are overlapping @ %d\n", index);

			i_begs++;
			score++;
		}

		//printf("Added begins together, score: %d\n", score);
	}
	else if((i_begs >= (int)frag_begs->size()) ||
		(frag_begs->at(i_begs) > frag_ends->at(i_ends)))
	{
		index = frag_ends->at(i_ends);
                i_ends++;
                score--;

		// Go over all the ends that have the same index in the ends list.
                while(i_ends < (int)frag_ends->size() &&
                        index == frag_ends->at(i_ends))
                {
			//printf("Ends are overlapping @ %d\n", index);
                        i_ends++;
                        score--;
                }

                //printf("Added ends together, score: %d\n", score);
	}
	else if((i_begs < (int)frag_begs->size()) &&
		(frag_begs->at(i_begs) == frag_ends->at(i_ends)))
	{
		// These two fragment ends cancel each other and the 
                index = frag_ends->at(i_ends);
		i_begs++;
		score++;
                while(i_begs < (int)frag_begs->size() &&
                        index == frag_begs->at(i_begs))
                {
                        i_begs++;
                        score++;
                }

		score--;
		i_ends++;
                while(i_ends < (int)frag_ends->size() &&
                        index == frag_ends->at(i_ends))
                {
                        i_ends++;
                        score--;
                }
	}

	// A sanity check: If the ends are consumed more than begins, there is a problem.
	if(i_begs < i_ends)
	{
		printf("The indices are messed up: %d (%d), %d (%d)\n", i_begs, frag_begs->at(i_begs), i_ends, frag_ends->at(i_ends));
		exit(0);
	}
}

/*
Count peaks for a list of thresholds. This function is more efficient than calling single threshold version since it processes the height 
data for all the thresholds in one go. Note also that the function counts the peaks within the region specified in the arguments.

This function is based on counting the peaks that are ended. Note that original PeakSeq implementation counts the peaks that are started,
but this implementation counts the peaks that overlap with a certain region, and needs to know exact start and end limits for the peaks
so as to be able to determine if there is an overlap between the requested region and the peak.
*/
void t_enrichment_profile::count_peaks(int* n_peaks_per_thresh, 
									   int min_thresh, 
									   int max_thresh, 
									   int min_gap_bw_peaks, 
									   int min_i_nuc, 
									   int max_i_nuc)
{
	signed int* next_possible_peak_starts = new signed int[max_thresh - min_thresh + 1];
	bool* above_threshold = new bool[max_thresh - min_thresh + 1];
	int* peak_starts = new int[max_thresh - min_thresh + 1];
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		next_possible_peak_starts[thresh - min_thresh] = -10000;
		n_peaks_per_thresh[thresh - min_thresh] = 0;
		above_threshold[thresh - min_thresh] = false;
		peak_starts[thresh - min_thresh ] = -10000;
	} // thresh loop.

    // Go over the profile and count the peaks.
    signed int height = 0;
	int cur_pos = 0; // The current fragment end position. This is the position for which the enrichment pfofile is last updated.
    for(int i_end = 0; i_end < (int)this->frag_ends->size(); i_end++)
    {
        // Update the current height: Note that there can be multiple ends at this position.
        height += this->frag_ends->at(i_end)->side;
        cur_pos = this->frag_ends->at(i_end)->pos;

        // Now add all the other ends at this position.
        while((i_end + 1) < (int)this->frag_ends->size() &&
                 cur_pos == this->frag_ends->at(i_end+1)->pos)
        {
            i_end++;
            height += this->frag_ends->at(i_end)->side;
        }

        if(height < 0)
        {
            fprintf(stderr, "Height below 0: %d (%s, %d)\n", height, __FILE__, __LINE__);
            exit(0);
        }

		// At this point, the height is changed. 
		// For each threhold that is lower than this height, this is either a new peak (if signal was below the threshold) otherwise 
		// must check if the closest peak was at least 200 bps away.
		for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
		{
			if(height >= thresh)
			{
				if(above_threshold[thresh - min_thresh])
				{
					// The signal was above the threshold in the last signal update, thus it has to be above this threshold. Do not need to do anything at this point.
				}
				else // Was below the threshold.
				{
					// check if this peak is far enough from the last peak.
					if(cur_pos > next_possible_peak_starts[thresh - min_thresh])
					{
						// This starts a new peak, and ends the previous pesk. Check if the last peak is in the requested window:
						int last_peak_end = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
						int last_peak_start = peak_starts[thresh - min_thresh];

						// If there is overlap between the peak and the investigated region, update the count by one.
						if(last_peak_start > max_i_nuc ||
							last_peak_end < min_i_nuc)
						{
						}
						else
						{
							n_peaks_per_thresh[thresh - min_thresh]++;
	                                                //printf("Added %d<%d-%d>\n", thresh, last_peak_start, last_peak_end);
						}

						// Mark the start of a new peak for this threshold.
						peak_starts[thresh - min_thresh] = cur_pos;
					}
				}
	
				// Set the new signal state as above threshold for this peak.
				above_threshold[thresh - min_thresh] = true;
			}
			else
			{			
				// If the signal was above this threshold previously, then this is the end of a peak, mark the next possible peak start for this trheshold.
				if(above_threshold[thresh - min_thresh])
				{
					next_possible_peak_starts[thresh - min_thresh] = cur_pos + min_gap_bw_peaks;
				}

				above_threshold[thresh - min_thresh] = false;
			}
		} // thresh loop
	} // i_end loop.

	// Ensure that the last peak is accounted for:
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		if(!above_threshold[thresh - min_thresh])
		{
			int last_peak_end = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
			int last_peak_start = peak_starts[thresh - min_thresh];

			if(last_peak_start > max_i_nuc ||
				last_peak_end < min_i_nuc)
			{
			}
			else
			{
				n_peaks_per_thresh[thresh - min_thresh]++;
			}
		}
		else
		{
			// The peak was above threshold in the last position. This means that all the positions remaining after this
			// will be incremented by one.
			if(peak_starts[thresh - min_thresh] > max_i_nuc)
			{
			}
			else
			{
				n_peaks_per_thresh[thresh - min_thresh]++;
			}
		}

                //printf("%d: %d\n", thresh, n_peaks_per_thresh[thresh - min_thresh]);
	} // i_thresh loop.

	// Free the array memory.
    delete [] next_possible_peak_starts;
    delete [] above_threshold;
    delete [] peak_starts;
}

t_enrichment_profile::~t_enrichment_profile()
{
	if(this->frag_ends == NULL)
	{
		// Dont do anything.
	}
	else if((int)this->frag_ends->size() == 0)
	{
		delete(this->frag_ends);
	}
	else
	{
		for(int i_frag_end = 0; i_frag_end < (int)this->frag_ends->size(); i_frag_end++)
		{
			free(this->frag_ends->at(i_frag_end));
		} // i_frag_end loop
	
		this->frag_ends->clear();
		delete(this->frag_ends);
	}
}

/*
int t_enrichment_profile::count_peaks(int threshold, int min_gap_bw_peaks, int min_i_nuc, int max_i_nuc)
{
        // Go over the profile and count the peaks.
        int peak_cnt = 0;
        signed int height = 0;
        signed int next_possible_peak_start = -10000;
	bool above_threshold = false;
        for(int i_end = 0; i_end < this->frag_ends->size(); i_end++)
        {
		// Update the current height: Note that there can be multiple ends at this position.
		height += this->frag_ends->at(i_end)->side;
		int cur_pos = this->frag_ends->at(i_end)->pos;

		// MUST CHECK FOR THE CASES WHEN THE PROFILE JUMPS OVER THE WHOLE REGION THAT IS BEING SEARCHED.

		// Now add all the other ends at this position.
		while((i_end + 1) < this->frag_ends->size() &&
			 cur_pos == this->frag_ends->at(i_end+1)->pos)
		{
			i_end++;
			height += this->frag_ends->at(i_end)->side;
		}

		if(height < 0)
		{
			fprintf(stderr, "Height below 0.\n");
			exit(0);
		}

		if(cur_pos >= min_i_nuc &&
			cur_pos < max_i_nuc)
		{
			// There are two conditions for this height to be a new peak:
			// 1. That it is higher than the threshold
			// 2. That its distance to last observed peak is larger than min_gap_bw_peaks.
	                if(height >= threshold)
        	        {
				if(above_threshold)
				{
				}
				else
				{
		                        if(cur_pos > next_possible_peak_start)
		                        {
						//printf("A peak @ %d\n", cur_pos);
		                                peak_cnt++;
		                        }
				}
				above_threshold = true;
	                }
	                else // Threshold is lower: There are two possibilities; 1. The previous region was a peak region, 2. It was not.
	                {
				if(above_threshold)
				{
	                                next_possible_peak_start = cur_pos + min_gap_bw_peaks;		
				}
				above_threshold = false;
	                }
		} // border check for nucleotides.
        } // i_end loop.

        return(peak_cnt);
}


// Count the peaks for this enrichment factor with the threshold.
int t_enrichment_profile::count_peaks(int threshold, int min_gap_bw_peaks)
{
	int min_i_nuc = 0;	
	int max_i_nuc = this->frag_ends->at(this->frag_ends->size() - 1)->pos + 1;
	return(this->count_peaks(threshold, min_gap_bw_peaks, min_i_nuc, max_i_nuc));
}
*/
/*
vector<t_peak_region*>* t_enrichment_profile::get_peaks(int threshold, int min_i_nuc, int max_i_nuc,  int min_gap_bw_peaks)
{
	//vector<t_peak_region*>* 
        // Go over the profile and count the peaks.
        int peak_cnt = 0;
        signed int height = 0;
        signed int next_possible_peak_start = -10000;
        bool above_threshold = false;
        for(int i_end = 0; i_end < this->frag_ends->size(); i_end++)
        {
                // Update the current height: Note that there can be multiple ends at this position.
                height += this->frag_ends->at(i_end)->side;
                int cur_pos = this->frag_ends->at(i_end)->pos;

                // MUST CHECK FOR THE CASES WHEN THE PROFILE JUMPS OVER THE WHOLE REGION THAT IS BEING SEARCHED.

                // Now add all the other ends at this position.
                while((i_end + 1) < this->frag_ends->size() &&
                         cur_pos == this->frag_ends->at(i_end+1)->pos)
                {
                        i_end++;
                        height += this->frag_ends->at(i_end)->side;
                }

                if(height < 0)
                {
                        fprintf(stderr, "Height below 0.\n");
                        exit(0);
                }

                if((cur_pos >= min_i_nuc &&
                        cur_pos < max_i_nuc) ||
			above_threshold)
                {
                        // There are two conditions for this height to be a new peak:
                        // 1. That it is higher than the threshold
                        // 2. That its distance to last observed peak is larger than min_gap_bw_peaks.
                        if(height >= threshold)
                        {
                                if(above_threshold)
                                {
                                }
                                else
                                {
                                        if(cur_pos > next_possible_peak_start)
                                        {
						//printf("A peak end @ %d>\n", next_possible_peak_start - min_gap_bw_peaks);
                                                //printf("<A peak start @ %d-", cur_pos);
                                                peak_cnt++;
                                        }
                                }
                                above_threshold = true;
                        }
                        else // Threshold is lower: There are two possibilities; 1. The previous region was a peak region, 2. It was not.
                        {
                                if(above_threshold)
                                {
                                        next_possible_peak_start = cur_pos + min_gap_bw_peaks;
                                }
                                above_threshold = false;
                        }
                } // border check for nucleotides.
        } // i_end loop.

	//printf("A peak end @ %d>\n", next_possible_peak_start - min_gap_bw_peaks);	

	return(NULL);
}
*/

/*
peaks_per_window_per_thresh[i_thresh][i_win] is the number of peaks in i_win'th for i_thresh.
*/

int** count_peaks_per_thresh_per_window(int n_meg_wins, 
					vector<t_peak_region*>** peaks_per_threshold,
					int min_thresh, 
					int max_thresh)
{
	// Allocate and initialize the peak counts per window per thresholds.
	int** peaks_per_window_per_thresh = new int*[max_thresh - min_thresh + 1];
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		peaks_per_window_per_thresh[thresh - min_thresh] = new int[n_meg_wins + 3];

		for(int i_win = 0; i_win < n_meg_wins; i_win++)
		{
			peaks_per_window_per_thresh[thresh - min_thresh][i_win] = 0;
		} // i_win loop 
	} // thresh loop.

	// Go over all the thresholds and all the peaks and update the counts for windows which overlaps with a peak.
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		// For this threshold, for each the peaks, update the count for the windows that overlap with the peak.
		vector<t_peak_region*>* cur_threshold_peaks = peaks_per_threshold[thresh - min_thresh];

		for(int i_peak = 0; i_peak < (int)cur_threshold_peaks->size(); i_peak++)
		{
			int start_i_win = cur_threshold_peaks->at(i_peak)->start_i_nuc / MEG_BASE;
			int end_i_win = cur_threshold_peaks->at(i_peak)->end_i_nuc / MEG_BASE;

			// update the counts for all windows that are overlapping with the peak.
			for(int i_win = start_i_win; i_win <= end_i_win; i_win++)
			{
				peaks_per_window_per_thresh[thresh - min_thresh][i_win]++;
			} // i_win loop
		} // i_peak loop.
	} // thresh loop.

	return(peaks_per_window_per_thresh);
}

void free_n_peaks_per_thresh_per_window(int** n_peaks_per_thresh_per_window, 
					int min_thresh,
                                        int max_thresh)
{
        //int** peaks_per_window_per_thresh = new int*[max_thresh - min_thresh + 1];
        for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
        {
                delete [] n_peaks_per_thresh_per_window[thresh - min_thresh];
        } // thresh loop.

	delete [] n_peaks_per_thresh_per_window;
}

void t_enrichment_profile::dump_peak_profiles_per_list(vector<t_peak_region*>* peak_region_list, char* dump_dir, char* file_name_fmt_str)
{
    // Go over the profile and count the peaks.
    signed int height = 0;
        int cur_pos = 0;
        int i_peak = 0;
        int i_end = 0;
    //for(int i_end = 0; i_end < this->frag_ends->size(); i_end++)
    while(i_end < (int)this->frag_ends->size() &&
                //i_peak < peak_region_list->size())
				i_peak < 1000)
    {
        //printf("@ peak %d <%d-%d>\n", i_peak, peak_region_list->at(i_peak)->start_i_nuc, peak_region_list->at(i_peak)->end_i_nuc);

        cur_pos = this->frag_ends->at(i_end)->pos;

        // Find the starting position for next peak and update the height along the way.
        while(cur_pos < peak_region_list->at(i_peak)->start_i_nuc)
        {
                // Update the current height: Note that there can be multiple ends at this position.
                height += this->frag_ends->at(i_end)->side;
                cur_pos = this->frag_ends->at(i_end)->pos;

                // Now add all the other ends at this position.
                while((i_end + 1) < (int)this->frag_ends->size() &&
                         cur_pos == this->frag_ends->at(i_end+1)->pos)
                {
                        i_end++;
                        height += this->frag_ends->at(i_end)->side;
                }

                i_end++;
        } // search for the starting of next peak.

	char cur_peak_prof_fp_fmt_str[1000];
	sprintf(cur_peak_prof_fp_fmt_str, "%s/%s", dump_dir, file_name_fmt_str);
	char cur_peak_prof_fp[1000];
	//sprintf(cur_peak_prof_fp, cur_peak_prof_fp_fmt_str, peak_region_list->at(i_peak)->start_i_nuc, peak_region_list->at(i_peak)->end_i_nuc);
	sprintf(cur_peak_prof_fp, cur_peak_prof_fp_fmt_str, i_peak);
	FILE* f_cur_peak_prof = open_f(cur_peak_prof_fp, "w");

	// If file is not open, return.
	if(f_cur_peak_prof == NULL)
	{
		printf("Could not open %s for writing peak profile data.\n", cur_peak_prof_fp);
		return;
	}

        fprintf(f_cur_peak_prof, "%d %d\n", cur_pos, height);


        while(cur_pos < peak_region_list->at(i_peak)->end_i_nuc)
        {
                // Update the current height: Note that there can be multiple ends at this position.
                height += this->frag_ends->at(i_end)->side;
                cur_pos = this->frag_ends->at(i_end)->pos;

                // Now add all the other ends at this position.
                while((i_end + 1) < (int)this->frag_ends->size() &&
                         cur_pos == this->frag_ends->at(i_end+1)->pos)
                {
                        i_end++;
                        height += this->frag_ends->at(i_end)->side;
                }

		fprintf(f_cur_peak_prof, "%d %d\n", cur_pos, height);

                i_end++;
        } // loop for processing current peak.

	fclose(f_cur_peak_prof);

	// Move to the next peak.
	i_peak++;
    }

}

// Set the maximum and the position of maximum for all the peak regions in the list using the enrichment profile.
void t_enrichment_profile::set_max_per_peak_region_list(vector<t_peak_region*>* peak_region_list)
{
    // Go over the profile and count the peaks.
    signed int height = 0;
        int cur_pos = 0;
	int i_peak = 0;
	int i_end = 0;
    //for(int i_end = 0; i_end < this->frag_ends->size(); i_end++)
    while(i_end < (int)this->frag_ends->size() &&
		i_peak < (int)peak_region_list->size())
    {
	//printf("@ peak %d <%d-%d>\n", i_peak, peak_region_list->at(i_peak)->start_i_nuc, peak_region_list->at(i_peak)->end_i_nuc);

        cur_pos = this->frag_ends->at(i_end)->pos;

	// Find the starting position for next peak and update the height along the way.
	while(cur_pos < peak_region_list->at(i_peak)->start_i_nuc)
	{
	        // Update the current height: Note that there can be multiple ends at this position.
	        height += this->frag_ends->at(i_end)->side;
	        cur_pos = this->frag_ends->at(i_end)->pos;
	
	        // Now add all the other ends at this position.
	        while((i_end + 1) < (int)this->frag_ends->size() &&
	                 cur_pos == this->frag_ends->at(i_end+1)->pos)
	        {
	                i_end++;
	                height += this->frag_ends->at(i_end)->side;
	        }
		
		i_end++;
	} // search for the starting of next peak.

	// Set the maximum and position for the maximum for this peak.
	bool height_at_max = true;
	peak_region_list->at(i_peak)->max_height = height;
	peak_region_list->at(i_peak)->max_height_start_i = cur_pos;

	while(cur_pos < peak_region_list->at(i_peak)->end_i_nuc)
	{
                // Update the current height: Note that there can be multiple ends at this position.
                height += this->frag_ends->at(i_end)->side;
                cur_pos = this->frag_ends->at(i_end)->pos;

                // Now add all the other ends at this position.
                while((i_end + 1) < (int)this->frag_ends->size() &&
                         cur_pos == (int)this->frag_ends->at(i_end+1)->pos)
                {
                        i_end++;
                        height += this->frag_ends->at(i_end)->side;
                }

		// Update maximum if necessary.
		if(peak_region_list->at(i_peak)->max_height < height)
		{
			peak_region_list->at(i_peak)->max_height = height;
			peak_region_list->at(i_peak)->max_height_start_i = cur_pos;
			height_at_max = true;
		}
		else if(height_at_max && peak_region_list->at(i_peak)->max_height > height)
		{
			peak_region_list->at(i_peak)->max_height_end_i = cur_pos - 1;
			height_at_max = false;
		}

		i_end++;
	}

	// Process the end of maximum if peak was still at maximum.
	if(height_at_max)
	{
		peak_region_list->at(i_peak)->max_height_end_i = peak_region_list->at(i_peak)->end_i_nuc;
		height_at_max = false;
	}

	// Move to the next peak.
	i_peak++;
    } // i_end loop.
}

/*
Determine all the peaks for this enrichment profile: There are two functions for this, one of them 
takes scaling_factor as a parameter to do thresholding on the scaled profile, that function uses
double values for profile. Following function uses integer values for storing the profile and does not
do scaling of the profile.
*/
vector<t_peak_region*>** t_enrichment_profile::get_peaks_per_thresholds(int min_thresh,
																		int max_thresh,
																		int min_gap_bw_peaks)
{
	// Allocate the peak region lists, a list per threshold is allocated.
	vector<t_peak_region*>** peak_regions_per_threshold = new vector<t_peak_region*>*[max_thresh - min_thresh + 1];
	for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
	{
		peak_regions_per_threshold[thresh - min_thresh] = new vector<t_peak_region*>();	
	} // thresh loop.	

    signed int* next_possible_peak_starts = new signed int[max_thresh - min_thresh + 1];
    bool* above_threshold = new bool[max_thresh - min_thresh + 1];
	bool* height_at_max = new bool[max_thresh - min_thresh + 1];
    for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
    {
            next_possible_peak_starts[thresh - min_thresh] = -10000;
            above_threshold[thresh - min_thresh] = false;
		height_at_max[thresh - min_thresh] = false;
    } // thresh loop.

    // Go over the profile and count the peaks.
    int height = 0;
	int cur_pos = 0;
    for(int i_end = 0; i_end < (int)this->frag_ends->size(); i_end++)
    {
        bool height_changed = false;
		int prev_height = height;

		// Update the current height: Note that there can be multiple ends at this position.
		height += this->frag_ends->at(i_end)->side;
		cur_pos = this->frag_ends->at(i_end)->pos;

		// Now add all the other ends at this position.
		while((i_end + 1) < (int)this->frag_ends->size() &&
					cur_pos == this->frag_ends->at(i_end+1)->pos)
		{
				i_end++;	
				height += this->frag_ends->at(i_end)->side;
		}

		// Check if the height is changed or not. If height is not changed, there is no need to process these fragment ends.
		if(prev_height == height)
		{
			height_changed = false;
		}
		else
		{
			height_changed = true;
		}

        if(height < 0)
        {
                fprintf(stderr, "Height below 0.\n");
                exit(0);
        }

		// At this point, the height is changed.
		// For each threhold that is lower than this height, this is either a new peak (if signal was below the threshold) otherwise
		// must check if the closest peak was at least 200 bps away.
		for(int thresh = min_thresh; 
			height_changed && 
			thresh <= max_thresh; 
			thresh++)
		{
			if(height >= thresh)
			{
				if(above_threshold[thresh - min_thresh])
				{
					// The signal was above the threshold in the last signal update, thus it has to be above this threshold. Do not need to do anything at this point.

					// Process the maximum height for this peak, if it is necessary.
					if((int)height > peak_regions_per_threshold[thresh - min_thresh]->back()->max_height)
					{
						peak_regions_per_threshold[thresh - min_thresh]->back()->max_height = (int)height;
						peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_start_i = cur_pos;
						height_at_max[thresh - min_thresh] = true;
					}

					// Check if the height fell below the maximum.
					if(height_at_max[thresh - min_thresh] && (height < peak_regions_per_threshold[thresh - min_thresh]->back()->max_height))
					{
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_end_i = cur_pos - 1;
                                                height_at_max[thresh - min_thresh] = false;			
					}
				}
				else // Was below the threshold.
				{
					// check if this peak is far enough from the last peak.
					if(cur_pos > next_possible_peak_starts[thresh - min_thresh])
					{
						//printf("Peak @ %d\n", cur_pos);

						// Set the ending of the previous peak.
						if((int)peak_regions_per_threshold[thresh - min_thresh]->size() > 0)
						{
							peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
							// Also if the peak was in the max height, set it as the end of the peak.
							if(height_at_max[thresh - min_thresh])
							{
								peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_end_i = cur_pos - 1;
								height_at_max[thresh - min_thresh] = false;
							}
						}

						t_peak_region* new_peak_region = new t_peak_region(cur_pos, 0); // The end is not initialized, yet.
						peak_regions_per_threshold[thresh - min_thresh]->push_back(new_peak_region); // Add this peak region to the list of peaks.
						// Set the maximum height and position of the height for this peak when the peak is initiated.
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height = (int)height;
                                                peak_regions_per_threshold[thresh - min_thresh]->back()->max_height_start_i = cur_pos;
                                                height_at_max[thresh - min_thresh] = true;
					}
				}

				// Set the new signal state as above threshold for this peak.
				above_threshold[thresh - min_thresh] = true;
			}
			else
			{
				// If the signal was above this threshold previously, then this is the end of a peak, mark the next possible peak start for this trheshold. Note that the maximum for the peak is not processed as the maximum height for a peak is always guaranteed to be at a position where the signal is higher than the threshold.
				if(above_threshold[thresh - min_thresh])
				{
					next_possible_peak_starts[thresh - min_thresh] = cur_pos + min_gap_bw_peaks-1;
				}
				above_threshold[thresh - min_thresh] = false;
			}
		} // thresh loop
	} // i_end loop.


	// Must check whether there is any peak that is hanging at the end. This is becaues a peak added when the next peak is starting.
	// Note that if there were no peaks in the whole chromosome, this loop should not be executed.
	for(int thresh = min_thresh; 
		thresh <= max_thresh;
		thresh++)
	{
		if((int)peak_regions_per_threshold[thresh - min_thresh]->size() > 0)
		{
			// If the last peak was below the threshold, this means that an end for this peak was found but the algorithm was looking for another peak to merge it with. Set the end for this peak.
			if(!above_threshold[thresh - min_thresh])
			{
				peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = next_possible_peak_starts[thresh - min_thresh] - min_gap_bw_peaks;
			}
			else
			{
				 peak_regions_per_threshold[thresh - min_thresh]->back()->end_i_nuc = cur_pos;				
			}
		}
	} // thresh loop

	return(peak_regions_per_threshold);	
}

void free_peaks_per_threshold(vector<t_peak_region*>** peak_regions_per_threshold,
				int min_thresh,
                              	int max_thresh)
{
        for(int thresh = min_thresh; thresh <= max_thresh; thresh++)
        {
                delete(peak_regions_per_threshold[thresh - min_thresh]);
        } // thresh loop.

	delete [] peak_regions_per_threshold;

}


