#ifndef __INDEXING_PER_FORMAT__
#define __INDEXING_PER_FORMAT__

/*
This file lists the indexing bases for different sequence formats. This becomes necessary when output of one method is used as input for another 
method. Note that some
*/

#define LIB_BASE		(1) // This is the indexing base for genomics libraries. All the libraries use 1 based indexing to refer to sequences.
#define ELAND_BASE		(1)
#define tagAlign_BASE	(1)
#define bowtie_BASE		(0) // Default bowtie output.
#define SAM_BASE		(0)
#define PEAKSEQ_BASE	(1)

#endif // __INDEXING_PER_FORMAT__

