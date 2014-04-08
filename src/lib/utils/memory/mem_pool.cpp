#include <stdio.h>
#include <stdlib.h>
#include "mem_pool.h"
#include <vector>
#include <string.h>

t_mem_pool::t_mem_pool()
{
	this->mem_pools = new vector<char*>();
	this->realloc_pool();
}

t_mem_pool::~t_mem_pool()
{
	// Free all the pools.
	for(int i_pool = 0; i_pool < this->mem_pools->size(); i_pool++)
	{
		delete[] this->mem_pools->at(i_pool);
	} // i_pool loop

	delete(this->mem_pools);
}

void* t_mem_pool::get_mem_from_pool(int size)
{
	// Is the usage going over the pool size?
	if(this->last_pool_mem_usage + size >= MEM_POOL_SIZE)
	{
		this->realloc_pool();
	}

	void* ptr = (void*)(this->current_pool + this->last_pool_mem_usage);
	this->last_pool_mem_usage += size;
	return(ptr);
}

void t_mem_pool::realloc_pool()
{
	// Allocate a new pool.
	printf("Re-allocing a new pool of %d bytes.\n", MEM_POOL_SIZE);
	this->last_pool_mem_usage = 0;
	this->current_pool = new char[MEM_POOL_SIZE];
	this->mem_pools->push_back(this->current_pool);

	// Initialize the memory to commit it.
	for(int i = 0; i < MEM_POOL_SIZE; i++)
	{
		this->current_pool[i] = 0;
	}
}


