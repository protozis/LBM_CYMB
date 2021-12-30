#include<stdio.h>
#include<stdlib.h>
#include"lbm.h"
#include"lbmcl.h"

void print_usage(){
	printf("Usage: laucher [config]\n");
	printf("Examples:\n");
	printf("\tlauncher a.conf\n");
}
int main(int argc, char *argv[]) {

	char confFileName[80];
	if(argc < 2){
		print_usage();
		fprintf(stderr,"launcher: missing configure file\n");
		exit(0);
	}
	sscanf(argv[1],"%s",confFileName);

	update_config(confFileName);
	simulate_ocl(ND_FILE,BC_FILE,OUTPUT_DIR,PROGRAM_FILE);
	return 0;
}


