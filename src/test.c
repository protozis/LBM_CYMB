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
	printf("testing\n");
	printf("open file %s\n",OBJ_FILE);
	FILE *f;
	struct OBJ *obj = NULL;
	struct CY *tmp = NULL;
	if(f = fopen(OBJ_FILE,"r")){
		obj = OBJ_read(f);
		fclose(f);
		printf("%lf\n",obj->m[0]->rist);
	}else{
		printf("can't open file\n");
	}

	return 0;
}


