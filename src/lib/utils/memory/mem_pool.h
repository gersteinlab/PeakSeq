#ifndef __MEM_POOL__
#define __MEM_POOL__

#include <vector>
using namespace std;

// The memory chunk size is 200 megs, corresponding to each memory allocation.
#define MEM_POOL_SIZE (200 * 1024 * 1024)

class t_mem_pool
{
public:
	t_mem_pool();
	~t_mem_pool();

	// The list of base pointers to the memory pools.
	vector<char*>* mem_pools;

	// This is the last allocated pool.
	char* current_pool;

	// This is the next pointer to return from the pool.
	int last_pool_mem_usage;

	// Return 
	void* get_mem_from_pool(int size);
	void realloc_pool();
};

#endif // __MEM_POOL__



