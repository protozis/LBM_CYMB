#include<stdio.h>
#include<stdlib.h>
#include"lbm.h"

void print_usage(){
	printf("Usage: nd_gen [OPTION]\n");
	printf("Otions:\n");
	printf("\t-x #, Numer of column\n");
	printf("\t-y #, Number of row\n");
	printf("\t-d #, Density (Default: 0.3)\n");
	printf("\t-i #, Ux (Default: 0.12)\n");
	printf("\t-j #, Uy (Default: 0)\n");
	printf("\t-h , Print this help page\n");

	printf("Example:\n");
	printf("\tnd_gen -x 100 -y 30 -i 0.06\n");
}



void main(int argc, char *argv[]) {

	int x = 0;
	int y = 0;
	double D =0.3;
	double ux=0.12;
	double uy=0;
	double *p;
	if (argc < 2){
		print_usage();
		exit(0);
	}
	while ((++argv)[0])
	{
		if (argv[0][0] == '-') {
			switch (argv[0][1]) {
				case 'x':
					sscanf(argv[1],"%d",&x);
					(++argv)[0];
					--argc;
					break;
				case 'y':
					sscanf(argv[1],"%d",&y);
					(++argv)[0];
					--argc;
					break;
				case 'd':
					sscanf(argv[1],"%lf",&D);
					(++argv)[0];
					--argc;
					break;
				case 'i':
					sscanf(argv[1],"%lf",&ux);
					(++argv)[0];
					--argc;
					break;
				case 'j':
					sscanf(argv[1],"%lf",&uy);
					(++argv)[0];
					--argc;
					break;
				case 'h':
					print_usage();
					exit(0);
					break;
			}
			--argc;

		}
		else{
			fprintf(stderr,"nd_gen: invalid argument -- '%s'\n",argv[0]);
			print_usage();
			exit(0);
			break;
		}
	}
	ND_init(p,D,ux,uy);
	printf("%d %d 9\n",x,y);
	for(int j=0;j<x*y;j++){
		for(int i=0;i<9;i++){
			printf("%lf ",p[i]);
		}
		printf("\n");
	}
}

