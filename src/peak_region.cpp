#include <stdio.h>
#include <stdlib.h>
#include "peak_region.h"

t_peak_region::t_peak_region()
{
	this->start_i_nuc = 0;
	this->end_i_nuc = 0;
	this->log_p_val = 0.0;
	this->q_val = 0.0;

	this->max_height = 0;
	this->max_height_start_i = 0;
    this->max_height_end_i = 0;
}

t_peak_region::t_peak_region(int _start_i_nuc, int _end_i_nuc)
{
	this->start_i_nuc = _start_i_nuc;
	this->end_i_nuc = _end_i_nuc;
        this->log_p_val = 0.0;
        this->q_val = 0.0;

        this->max_height = 0;
        this->max_height_start_i = 0;
        this->max_height_end_i = 0;
}

t_peak_region::t_peak_region(int _start_i_nuc, int _end_i_nuc, double _log_p_val)
{
        this->start_i_nuc = _start_i_nuc;
        this->end_i_nuc = _end_i_nuc;
        this->log_p_val = _log_p_val;
        this->q_val = 0.0;

        this->max_height = 0;
        this->max_height_start_i = 0;
        this->max_height_end_i = 0;
}

t_peak_region::~t_peak_region()
{
}

