#ifndef __CONFIG__
#define __CONFIG__

#include <vector>
using namespace std;

typedef vector<char*> t_val_list;

class t_config
{
public:
	t_config(char* config_fp);
	~t_config();

	// ID's and values for the entries in the configuration file.
	vector<char*>* ids;
	vector<t_val_list*>* vals; // A value list for each identifier.

	// Can oly access the entries using th id's.
	bool get_double_val(char* id, double& d_val);
	bool get_str_val(char* id, char* val_buff);
	bool get_int_val(char* id, int& int_val);
	t_val_list* get_val_list(char* id);
};

#endif // __CONFIG__
