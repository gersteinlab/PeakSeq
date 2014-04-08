PeakSeq1.1 
===========
PeakSeq is a program for calling peaks in ChIP-Seq experiments. Requires a set of mapped reads from ChIP enriched and Input DNA (or some other control) based
sequencing experiment.

Refer to the paper for explanation of the method
PeakSeq Systematic Scoring of ChIP-Seq Experiments Relative to Controls
Rozowsky J, Euskirchen G, Auerbach R, Zhang Z, Gibson T, Bjornson R, Carriero N, Snyder M, Gerstein M
Nature Biotechnology 27, 66 - 75 (2009).

This new version is re-coded by Arif Harmanci under supervision of Joel Rozowsky.

Installation 
=============
Extract the folder, then do make. PeakSeq is built into bin directory.

Usage
======
The mapped read files that contain the ChIP-Seq reads and the control reads must be preprocessed into a binary file format that PeakSeq can access fast. This is performed by running PeakSeq with -preprocess option. Currently PeakSeq can preprocess three mapped read file formats SAM, tagAlign, and ELAND. 

For example, following command reads the mapped reads in the file named pol2_chip_seq.eland and dumps the processed reads into the folder Pol2_ChIP
.PeakSeq -preprocess ELAND chr_id_list.txt pol2_chip_seq.eland Pol2_ChIP
chr_id_list.txt is the list of chromosomes to be processed. These must match chromosome ids in the mapped read file. 

BAM files can be processed first by samtools package then be piped into PeakSeq
samtools view pol2_chip_seq.bam  .PeakSeq -preprocess SAM chr_id_list.txt stdin Pol2_ChIP
where the input file name is specified as stdin to indicate that the mapped reads are piped into PeakSeq.

After preprocessing, the peaks can be called using -peak_select option. PeakSeq requires a configuration file for calling the peaks where several parameters are specified in each line in the '[id] [value]' format. An example configuration file is as follows

config.dat
# Experiment id is used as a prefix to the output file name.
Experiment_id pol2_peaks

# Chromosome ID list file, is used to generate data file paths.
chromosome_list_file chr_id_list.txt

# Enrichment fragment length For tag extension, this is the value of average fragment length.
Enrichment_fragment_length 200

# Target FDR in the simulations.
target_FDR 0.05

# Number of simulations performed while estimating the putative peaks.
N_Simulations 10

# Minimum distance between consecutive peaks
Minimum_interpeak_distance 200

# Mappability file that includes the uniquely mappable number of nucleotides per window for each chromosome.
Mappability_map_file Mapability_HG.txt

# The directory that contains the preprocessed ChIP-Seq reads, can specify multiple directories to pool reads from multiple source (e.g. replicates)
ChIP_Seq_reads_data_dirs Pol2_ChIP

# The directory that contains the preprocessed Input (control) experiment reads. (Multiple directories allowed)
Input_reads_data_dirs Pol2_Input

# Seed for pseudo-random number generator. This is necessary for simulated background option (specified below).
#Simulation_seed 1234567

# Q-value threshold applied on the final set of peaks.
max_Qvalue 0.05

# There are currently two models for simulating the background for threshold selection
# Simulated background is the simulation based method that is explained in the PeakSeq paper.
# Poisson background uses a simple Poisson background with mean estimated from the read statistics. This option is still experimental but it is much faster than the simulated background option.
# Background_model Poisson
Background_model Simulated


config.dat is passed to PeakSeq as a parameter with -peak_select option
PeakSeq -peak_select config.dat.

PeakSeq will complain if it cannot find a mandatory entry in the configuration file.

Output
=======
PeakSeq dumps a narrowPeak file (httpencodewiki.ucsc.eduEncodeDCCindex.phpFile_Formats#narrowPeak_Narrow_.28or_Point-Source.29_Peaks_Format) and a txt file with peaks. The peak region, p-value, q-value, and enrichment statistics are reported in the output file.

FeedbackCommentsQuestions
============================
Are always welcome. Please address them to Arif Harmanci (arif.harmanci@yale.edu). 

=========================================================================
Arif Harmanci, Joel Rozowsky. Gerstein Lab, Yale University. June, 2011.
=========================================================================

