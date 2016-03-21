#include <stdio.h>
#include <stdlib.h>
#include "rng.h"

#include <math.h>
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

t_rng::t_rng(int _prn_gen_seed)
{
	this->iff = 0; // not be modifed; see Knuth.
	this->prn_gen_seed = _prn_gen_seed;
}

t_rng::~t_rng()
{};

double t_rng::random_double_ran3()
{
	long idum = this->prn_gen_seed;
	double random_num = (double)this->ran3(&idum);
	//printf("Generated %f\n",random_num);
	return(random_num);
}

// idum number is the seed that is used to "warm up" the generator by 
// setting ma array. idum is set to 1 after that, i.e., it is not used any more.
// ma array is important and should not be changed while generating a sequence of random numbers.
// After initialization of ma array is utilized to generate random numbers.
double t_rng::ran3(long* idum)
{
	long mj,mk;
	int i,ii,k;

	// Initialization.
	if (*idum < 0 || iff == 0) 
	{ 
		iff=1;
		mj=labs(MSEED-labs(*idum)); // Initialize ma[55] using the seed idum and the
		mj %= MBIG; // large number MSEED.
		ma[55]=mj;
		mk=1;

		// Now initialize the rest of the table,
		for (i=1;i<=54;i++) 
		{ 
			ii=(21*i) % 55; // in a slightly random order,
			ma[ii]=mk; // with numbers that are not especially random.
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}

		// We randomize them by warming up the generator
		for (k=1;k<=4;k++) 
		{
			for(i=1;i<=55;i++) 
			{
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) 
				{
					ma[i] += MBIG;
				}
			}
		}

		inext=0; 
		inextp=31; // The constant 31 is special; see Knuth.
		*idum=1;
	}

	if (++inext == 56) 
	{
		inext=1; // Increment inext and inextp, wrapping around
	}

	if (++inextp == 56) 
	{
		inextp=1; 
	}

	// Generate a new random number subtractively.
	mj=ma[inext]-ma[inextp]; 

	// Be sure that it is in range.
	if (mj < MZ) 
	{
		mj += MBIG; 
	}

	ma[inext]=mj; // Store it,
	return mj*FAC; // and output the derived uniform deviate.
}


