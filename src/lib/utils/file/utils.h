#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>
#include <vector>

using namespace std;

#define __MAX_PATH (10000)

#define DATAPATH_ENV_VAR "DATAPATH"
#define LOCAL_DATA_PATH "data"

char* resolve_data_dir();
bool check_file(char* fp);
void validate_file(char* fp);
char* x_fgets(char* buff, int size, FILE* file);
FILE* open_f(char* fp, char* mode);
vector<char*>* load_directory_files(char* root_dir, char* extension);
char* get_file_name(char* fp);

vector<char*>* buffer_file(char* fp);

// Loads a line from a file, without the buffer size. In principle, reads line of any length and buffers it.
char* getline(FILE* file);

bool line_empty(char* line);

#endif // _UTILS_

