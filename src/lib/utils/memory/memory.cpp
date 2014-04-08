#include <stdio.h>
#include <stdlib.h>

#include "memory.h"
#include "../exception_obj/exception_obj.h"

/*
Exceptional malloc.
*/
void* alloc_mem(size_t size)
{	
	unsigned char* mem = NULL;
	mem = (unsigned char*)malloc(size);
	if(mem == NULL)
	{
		throw(new t_exception_obj("Memory allocation failed."));
	}

	return((void*)mem);
}
