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
	printf("open file %s\n",BC_FILE);
	FILE *f;
	struct BC *bc = NULL;
	if(f = fopen(BC_FILE,"r")){
		bc = BC_read(f);
		printf("%lf %lf\n",((struct CY *)bc->m[0])->pos[0],((struct CY *)bc->m[0])->pos[1]);
		printf("%lf %lf %lf\n",bc->dnt,bc->ux,bc->uy);
		fclose(f);
	}else{
		printf("can't open file\n");
	}

	return 0;
}


