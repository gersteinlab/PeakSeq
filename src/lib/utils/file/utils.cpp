#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

vector<char*>* load_directory_files(char* root_dir, char* extension)
{
        char ls_cmd[1000];

		if(extension != NULL)
		{
			sprintf(ls_cmd, "ls %s/*.%s > cmd_op.txt", root_dir, extension);
			system(ls_cmd);
		}
		else
		{
			sprintf(ls_cmd, "ls %s > cmd_op.txt", root_dir);
			system(ls_cmd);
		}

        vector<char*>* dir_files = new vector<char*>();

        FILE* f_cmd_op = fopen("cmd_op.txt", "r");
        char current_fn[1000];
        while(fscanf(f_cmd_op, "%s", current_fn) == 1)
        {
                char* new_fn = new char[strlen(current_fn) + 3];
                strcpy(new_fn, current_fn);
                dir_files->push_back(new_fn);
        }
        fclose(f_cmd_op);

		// Erase the file.
        sprintf(ls_cmd, "rm -f cmd_op.txt", root_dir);
        system(ls_cmd);

        return(dir_files);
}

vector<char*>* buffer_file(char* fp)
{
	vector<char*>* file_lines = new vector<char*>();
	printf("Buffering %s.\n", fp);

	int n_lines = 0;
	FILE* f = fopen(fp, "r");	
	while(1)
	{
		char* new_line = getline(f);
		if(new_line == NULL)
		{
			break;
		}
		file_lines->push_back(new_line);

		n_lines++;
		if((n_lines % 10000) == 0)
		{
			printf("Loaded %d lines.               \r");
		}
	} // file buffering loop.

	fclose(f);

	return(file_lines);
}

bool check_file(char* fp)
{
	FILE* f_temp = fopen(fp, "r");
	if(f_temp == NULL)
	{		
		return(false);
	}

	fclose(f_temp);
	return(true);
}

// Check for valid CR LF's depending on OS,
// Should be run for all ASCII input files.
#define CR (0x0D)
#define LF (0x0A)
void validate_file(char* fp)
{
#ifdef _WIN32
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		if(cur_char == CR)
		{
			if(fread(&cur_char, 1, 1, f_ip_bin) == 1)
			{
				if(cur_char != LF)
				{
					// Just a warning here.
					printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
					//exit(0);
				}
			}
			else
			{
				// Just a warning here.
				printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
				//exit(0);
			}
		}
		else if(cur_char == LF) // If there is an immediate LF before seeing a CR, this is a linux file.
		{
			// Just a warning here.
			printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}

	}
	fclose(f_ip_bin);
#endif

#ifdef __unix__
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// Linux files do not contain CR's.
		// They only contain LF's.
		if(cur_char == CR)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif
}

FILE* open_f(char* fp, char* mode)
{
	if(fp == NULL || mode == NULL)
	{
		printf("Invalid arguments to open_f.\n", fp);
		exit(0);
	}

	FILE* f = fopen(fp, mode);

	if(f == NULL)
	{
		if(mode[0] == 'r')
		{
			printf("Could not open %s for reading.\n", fp);
			exit(0);
		}
		else if(mode[0] == 'w')
		{
			printf("Could not open %s for writing.\n", fp);
			exit(0);
		}
		else
		{
			printf("Could not open %s for requested operation.\n", fp);
			exit(0);
		}
	}

	return(f);
}

bool line_empty(char* line)
{
	int l = strlen(line);
	for(int i = 0; i < l; i++)
	{
		if(line[i] != ' ' &&
			line[i] != '\t' &&
			line[i] != '\n')
		{
			return(false);
		}
	}

	return(true);
}

char* getline(FILE* file)
{
	vector<char>* line_vec = new vector<char>();

	char ret = 0;
	while(1)
	{
		ret = getc(file);

		if(ret == EOF)
		{
			break;
		}
		else if(ret == '\n')
		{
			break;
		}

		// Add this character to the vector of characters.
		line_vec->push_back(ret);
	} // file reading loop.

	// If the end-of-file is encountered and there was nothing read, return NULL to indicate the there is nothing left to read and EOF is reached.
	if(ret == EOF && 
		line_vec->size() == 0)
	{
		delete(line_vec);
		return(NULL);
	}

	char* buff = new char[line_vec->size() + 2];
	memset(buff, 0, sizeof(char) * (line_vec->size() + 2));

	// Copy the vector to buffer.
	for(int i = 0; i < line_vec->size(); i++)
	{
		buff[i] = line_vec->at(i);
	} // i loop.

	delete(line_vec);

	return(buff);
}

char* resolve_data_dir()
{
	// try to resolve the DATAPATH_ENV_VAR.
	char* data_dir_from_env = getenv(DATAPATH_ENV_VAR);

	if(data_dir_from_env != NULL)
	{
		return(data_dir_from_env);
	}
	else
	{
		return(LOCAL_DATA_PATH);
	}

	printf("Could not resolve thermodynamics data directory.\n");
	exit(0);
}

char* x_fgets(char* buff, int size, FILE* file)
{
	if(fgets(buff, size, file) == NULL)
	{
		return(NULL);
	}

	if(buff[strlen(buff) - 1] == '\n')
	{
		int i_new_line_char = (int)strlen(buff) - 1;
		buff[i_new_line_char] = 0;
	}

	return(buff);
}




