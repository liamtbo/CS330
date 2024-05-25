#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <linux/limits.h>

void usage(int argc, char** argv);
void find_file(char* dir_name, char* file_to_find);

int main(int argc, char** argv)
{
	DIR* dp;
	struct dirent* dirp;

	// check if this application is being used properly
	usage(argc, argv);

	// check to see if provided directory is accessible
	errno = 0;
	dp = opendir(argv[1]);
	if(dp == NULL) {
		switch(errno) {
			case EACCES:
				fprintf(stderr, "Permission denied\n");
				break;
			case ENOENT:
				fprintf(stderr, "Directory does not exist\n");
				break;
			case ENOTDIR:
				fprintf(stderr, "'%s' is not a directory\n", argv[1]);
				break;	
		}
	}

	// print all files in the directory
	int cnt = 0;
	while((dirp = readdir(dp)) != NULL) {
		fprintf(stdout, "%d: %s", cnt, dirp->d_name);
		if(dirp->d_type == DT_DIR) {
			printf("\t directory");
		}
		printf("\n");
		cnt++;
	}

	// close the directory 
	closedir(dp);


	// now, recursivey traverse the directory structure to find the provided
	// file name
	char* file_to_find = argv[2];
	find_file(argv[1], file_to_find);

	return 0;
}


void usage(int argc, char** argv)
{
    if (argc != 3) {
        fprintf(stderr, "Usage: ./%s directory_name file_to_find\n", argv[0]);
        exit(EXIT_FAILURE);
    }
}

void find_file(char* dir_name, char* file_to_find)
{
    // initializes
    DIR* dp; 
    struct dirent* dirp;
	dp = opendir(dir_name); // open directory

    // current directory curr_path
    char curr_path[256] = "";
    strcat(curr_path, dir_name);
    strcat(curr_path, "/");
    // printf("%s\n", curr_path);

    // loop over current dir contents
	while((dirp = readdir(dp)) != NULL) { // loopover every element in the dir
        
        // tmp path builder, updated in when dir element is a dict
        char path_builder[256] = "";
        strcat(path_builder, curr_path);
        // if directory, recursively search it
        if (dirp->d_type == DT_DIR) {
            // exclude directory . and .. in every sub dir
            if (strcmp(dirp->d_name, ".") != 0 && strcmp(dirp->d_name, "..") != 0) {
                
                // char *curr_path_builder = strcat(strcat(dir_name, "/"), dirp->d_name);
                strcat(path_builder, dirp->d_name);
                find_file(path_builder, file_to_find); 
            }
        } else { // if file, check if it's file_to_find
            if (strcmp(dirp->d_name, file_to_find) == 0) {
                printf("Found %s in %s\n", file_to_find, dir_name);
            }
        }
	}
    closedir(dp);
}