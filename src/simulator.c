#include<stdio.h>
#include<stdlib.h>
#include"lbm.h"
#include"lbmcl.h"

void print_usage(){
	printf("Usage: simulate_ocl [OPTION]... [ND_FILE]\n");
	printf("Options:\n");
	printf("\t-i [ND_FILE],\tNumber Density File\n");
	printf("\t-b [BC_FILE],\tBoundary Condition File\n");
	printf("\t-n, #,\tLOOP for #(int) times. (default: 1)\n");
	printf("\t-f, #,\tskip fram for #(int) times. (default: 1)\n");
	printf("\t-c #,\tcollision frequency(double) (default: 1)\n");
	printf("\t-s #,\tspeed of lattice(double) (default: 1)\n");
	printf("\t-p [PD_FILE],\tPlatform/device profile\n");
	printf("\t-o [DIR],\toutput directory\n");
	printf("\t-m,\tstream ND output(default: false)\n");
	printf("\t-u,\tND output to file(default: false)\n");
	printf("\t-q,\tSave final data\n");
	printf("\t-h,\tprint this help page\n");

	printf("Examples:\n");
	printf("\tsimulate_ocl -n 100 -i test.nd\n");
}




int main(int argc, char *argv[]) {

	if(argc < 2){
		print_usage();
		fprintf(stderr,"simulate_ocl: missing argument\n");
		exit(0);
	}
	while ((++argv)[0])
	{
		if (argv[0][0] == '-') {
			switch (argv[0][1]) {
				case 'n':
					sscanf(argv[1],"%d",&LOOP);
					(++argv)[0];
					--argc;
					break;
				case 'f':
					sscanf(argv[1],"%d",&SKP);
					(++argv)[0];
					--argc;
					break;
				case 'c':
					sscanf(argv[1],"%lf",&CF);
					(++argv)[0];
					--argc;
					break;
				case 's':
					sscanf(argv[1],"%lf",&SL);
					(++argv)[0];
					--argc;
					break;
				case 'o':
					sscanf(argv[1],"%s",DIR_NAME);
					(++argv)[0];
					--argc;
					break;
				case 'b':
					sscanf(argv[1],"%s",BC_FILE);
					(++argv)[0];
					--argc;
					break;
				case 'i':
					sscanf(argv[1],"%s",ND_FILE);
					(++argv)[0];
					--argc;
					break;
				case 'p':
					sscanf(argv[1],"%s",PD_FILE);
					(++argv)[0];
					--argc;
					break;
				case 'm':
					IS_MP4 = 1;
					break;
				case 'u':
					IS_FILE_OUTPUT = 1;
					break;
				case 'q':
					IS_SAVE_DATA = 1;
					break;
				case 'h':
					print_usage();
					exit(0);
					break;
			}
			--argc;

		}
		else{
			fprintf(stderr,"simulate_ocl: invalid argument -- '%s'\n",argv[0]);
			print_usage();
			exit(0);
			break;
		}
	}
	if(ND_FILE[0] == '\0'){
		fprintf(stderr,"simulate_ocl: missing ND file\n");
		exit(0);
	}
	if(BC_FILE[0] == '\0'){
		fprintf(stderr,"simulate_ocl: missing BC file\n");
		exit(0);
	}
	if(DIR_NAME[0] == '\0'){
		fprintf(stderr,"simulate_ocl: missing output directory\n");
		exit(0);
	}

	simulate_ocl(ND_FILE,BC_FILE,PD_FILE,DIR_NAME);

	return 0;
}


