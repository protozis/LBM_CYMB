#define KERNEL_FUNC "propagate"
#define CL_TARGET_OPENCL_VERSION 300
#define FC_OFFSET 1000
#define CS_LTTC 0.57735
#define PPMAX 255
#define PI 3.1415926


#define MP4_NUM 3

#ifdef MAC
#include<OpenCL/cl.h>
#else
#include<CL/cl.h>
#endif

#include"lbm.h"
#include"lbmcl.h"
#include<stdio.h>
#include<string.h>
#include<math.h>

//Default
size_t LOOP = 1;
size_t SKP = 1;
double CF = 1;
double CS = 340;
double CL = 1;
double CT = 1;
double CD = 1;
double MA = 0.1;
int IS_MP4 = 0;
int IS_SAVE_DATA = 0;
int IS_FILE_OUTPUT = 0;
char ND_FILE[80];
char BC_FILE[80];
char PD_FILE[80];
char OUTPUT_DIR[80];
char PROGRAM_FILE[80];
double REFUEL_RTO = 0.8;
double EAT_RTO = 0.01;
char LOG_FILE[80];
char DEBUG_FILE[80];
int IS_LOG_PRINT = 1;
int IS_PROGRESS_PRINT = 1;
double PL_MAX_D = 0.5;
double PL_MAX_UX = 0.1;
double PL_MAX_UY = 0.1;

double Cspring;
double Cdamp;
double Cmass;
double Cforce;
double Cacc;
double Cvel;
double Csl;


void update_config(char* filename){
	FILE *file;
	char par[80];
	char val[80];
	char *ptr;
	if(file = fopen(filename,"r")){
		while(fscanf(file,"%s %s\n",&par,&val) != EOF){
			if(strcmp(par,"LOOP") == 0){
				LOOP = atoi(val);
			}else if(strcmp(par,"SKP") == 0){
				SKP = atoi(val);
			}else if(strcmp(par,"CF") == 0){
				CF = strtod(val,&ptr);
			}else if(strcmp(par,"CS") == 0){
				CS = strtod(val,&ptr);
			}else if(strcmp(par,"CL") == 0){
				CL = strtod(val,&ptr);
			/*
			}else if(strcmp(par,"CT") == 0){
				CT = strtod(val,&ptr);
			*/
			}else if(strcmp(par,"CD") == 0){
				CD = strtod(val,&ptr);
			}else if(strcmp(par,"MA") == 0){
				MA = strtod(val,&ptr);
			}else if(strcmp(par,"IS_MP4") == 0){
				IS_MP4 = atoi(val);
			}else if(strcmp(par,"IS_SAVE_DATA") == 0){
				IS_SAVE_DATA = atoi(val);
			}else if(strcmp(par,"IS_FILE_OUTPUT") == 0){
				IS_FILE_OUTPUT = atoi(val);
			}else if(strcmp(par,"ND_FILE") == 0){
				 memcpy(ND_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"BC_FILE") == 0){
				 memcpy(BC_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"PD_FILE") == 0){
				 memcpy(PD_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"OUTPUT_DIR") == 0){
				 memcpy(OUTPUT_DIR,val,strlen(val)+1);
			}else if(strcmp(par,"PROGRAM_FILE") == 0){
				 memcpy(PROGRAM_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"REFUEL_RTO") == 0){
				REFUEL_RTO = strtod(val,&ptr);
			}else if(strcmp(par,"EAT_RTO") == 0){
				EAT_RTO = strtod(val,&ptr);
			}else if(strcmp(par,"LOG_FILE") == 0){
				 memcpy(LOG_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"DEBUG_FILE") == 0){
				 memcpy(DEBUG_FILE,val,strlen(val)+1);
			}else if(strcmp(par,"IS_LOG_PRINT") == 0){
				IS_LOG_PRINT= atoi(val);
			}else if(strcmp(par,"IS_PROGRESS_PRINT") == 0){
				IS_PROGRESS_PRINT = atoi(val);
			}else if(strcmp(par,"PL_MAX_D") == 0){
				PL_MAX_D = strtod(val,&ptr);
			}else if(strcmp(par,"PL_MAX_UX") == 0){
				PL_MAX_UX = strtod(val,&ptr);
			}else if(strcmp(par,"PL_MAX_UY") == 0){
				PL_MAX_UY = strtod(val,&ptr);
			}
		}

	}else{
		fprintf(stderr,"read_config: can't open conf file\n");
	}
}
void set_parameters(struct BC *bc,struct ND *nd){
	if(sqrt(pow(bc->ux,2) + pow(bc->uy,2)) == 0 ){
		Csl = 0.1/MA;
	}else{
		Csl = sqrt(pow(bc->ux,2) + pow(bc->uy,2)) / MA;
	}
	CT = CL*Csl/CS;
	Cmass = CD/(pow(CL,3));
	Cspring = Cmass/pow(CT,2);
	Cdamp = Cmass/CT;
	Cacc = CL/pow(CT,2);
	Cforce = Cmass*Cacc;
	Cvel = CL/CT;
}
void print_parameters(FILE *f,struct BC *bc,struct ND *nd){
	fprintf(f,"[Parameters]: (* config value)\n");
	fprintf(f,"\t<Stratage>\n");
	fprintf(f,"\t\t1. Similarity for the Reynolds number\n");
	fprintf(f,"\t\t2. Spectify CL for CT\n");
	
	fprintf(f,"\t<Dimensionless>\n");
	fprintf(f,"\t *\tMach number(MA): %lf\n",MA);
	fprintf(f,"\t\tReynolds number: %lf\n",nd->ny*MA/(CS_LTTC*((1/CF)-0.5)));
	fprintf(f,"\t\tGrid Reynolds number: %lf\n",bc->ux/(CS_LTTC*((1/CF)-0.5)));
	
	fprintf(f,"\t<Lattice unit>\n");
	fprintf(f,"\t *\tCollision frequency(CF): %lf\n",CF);
	fprintf(f,"\t\tKinematic viscosity: %lf\n",(1/CF)-0.5)/3;
	fprintf(f,"\t *\tBCV D: %lf Ux: %lf Uy: %lf\n",bc->dnt,bc->ux,bc->uy);
	fprintf(f,"\t *\tSize nx: %d ny: %d\n",nd->nx,nd->ny);
	fprintf(f,"\t\tSpeed of sound(Csl): %lf\n",Csl);

	fprintf(f,"\t<Conversion factors>\n");
	fprintf(f,"\t *\tLength(CL): %lf (m/lattice space)\n",CL);
	fprintf(f,"\t\tTime(CT): %lf (secs/time step)\n",CL*CS_LTTC/CS);
	fprintf(f,"\t *\tDensity(CD): %lfkg/m^3\n",CD);
	fprintf(f,"\t\tMass: %lfkg\n",Cmass);
	fprintf(f,"\t\tForce: %lfkg*m/s^2\n",Cforce);
	fprintf(f,"\t\tSpring constant: %lfkg/s^2\n",Cspring);
	fprintf(f,"\t\tDamping constant: %lfkg/s\n",Cdamp);

	fprintf(f,"\t<SI unit>\n");
	fprintf(f,"\t\tKinematic viscosity: %lfm^2/s\n",(((1/CF)-0.5)/3)*CL*CL/CT);
	fprintf(f,"\t\tSize width: %lfm height: %lfm\n",nd->nx*CL,nd->ny*CL);
	fprintf(f,"\t *\tSpeed of sound(CS): %lfm/s\n",CS);
	fprintf(f,"\t\tBCV D: %lfkg/m^3 Ux: %lfm/s Uy: %lfm/s\n",bc->dnt*CD,bc->ux*CL/CT,bc->uy*CL/CS);

	fprintf(f,"\t<Dirty tricks>\n");
	fprintf(f,"\t *\tREFUEL_RTO: %lf\n",REFUEL_RTO);
	fprintf(f,"\t *\tEAT_RTO: %lf\n",EAT_RTO);
	
	fprintf(f,"\t<Objects>\n");
	fprintf(f,"\t\t[spring] [damping] [mass] [Nau_freq] [Nau_cyc]\n");
	struct CY *tmp = NULL;
	double nf; 
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		nf = sqrt(tmp->spring/tmp->mass)/(2*PI*CT);
		fprintf(f,"\t\t%d: %lfkg/s^2 %lfkg/s %lfkg %lfHz %lfs\n",i,tmp->spring*Cspring,tmp->damp*Cdamp,tmp->mass*Cmass,nf,1/nf);
	}
}

void jet_colormap(double val, double min, double max,int *c){
	double dv = max - min;
	double dc[3];
	double sv;
	if(val < min){
		sv = -1;
	} else if (val > max){
		sv = 1;
	} else {
		sv = -1 + 2*(val - min)/dv;
	}

	if(sv < -0.75){
		dc[0] = 0;
		dc[1] = 0;
		dc[2] = 2*sv +2.5;
	} else if (sv < -0.25) {
		dc[0] = 0;
		dc[1] = 2*sv + 1.5;
		dc[2] = 1;
	} else if (sv < 0.25) {
		dc[0] = 2*sv + 0.5;
		dc[1] = 1;
		dc[2] = -2*sv + 0.5;
	} else if (sv < 0.75) {
		dc[0] = 1;
		dc[1] = -2*sv + 1.5;
		dc[2] = 0;
	} else {
		dc[0] = -2*sv + 2.5;
		dc[1] = 0;
		dc[2] = 0;
	}
	for(int i=0;i<3;i++){
		c[i] = (int)(PPMAX*dc[i]);
	}
}

void nd_ppm_write(double *m, int nx, int ny, double pmax, double pmin, int scale, FILE *f){
	fprintf(f,"P3\n");
	fprintf(f,"%d %d 255\n",nx*scale,ny*scale);
	int *v;
	for(int j=0;j<ny;j++){
		for(int sy=0;sy<scale;sy++){
			for(int i=0;i<nx;i++){
				for(int sx=0;sx<scale;sx++){
					v=(int *)malloc(3*sizeof(int));
					jet_colormap(m[i+j*nx],pmin,pmax,v);
					fprintf(f,"%d %d %d\n",v[0],v[1],v[2]);
					free(v);
				}
				int nx,ny;
			}
		}
	}
}
int check_exist(char *fn){
	FILE *f;
	if(f=fopen(fn,"r")){
		fclose(f);
		return 1;
	}else{
		return 0;
	}
}
void simulate_ocl(char* ndFileName, char* bcFileName, char* pdFileName, char* dirName, char* programFileName) {

	cl_device_id device;
	cl_context context;
	cl_program program;
	cl_kernel kernel;
	cl_command_queue queue;
	cl_int err;
	cl_mem nd_buffer;
	cl_mem res_buffer;
	cl_mem bcv_buffer;
	cl_mem bcpos_buffer;
	cl_mem bcpos_p_buffer;
	cl_mem bcvel_buffer;
	cl_mem bcrad_buffer;
	cl_mem bcfc_buffer;

	struct BC *bc = NULL;
	struct ND *nd = NULL;
	FILE *output;
	FILE *debug_out;
	FILE *log_out;
	char *log_line;
	size_t log_len;
	if(!(debug_out = fopen(DEBUG_FILE,"w"))){
		fprintf(stderr,"simulate_ocl: can't open debug file %s\n",DEBUG_FILE);
		exit(1);
	}
	if(!(log_out = fopen(LOG_FILE,"w"))){
		fprintf(stderr,"simulate_ocl: can't open log file %s\n",LOG_FILE);
		exit(1);
	}
	FILE *input;
	char filename[80];
	if(input = fopen(bcFileName,"r")){
		bc = BC_read(input);
	} else {
		fprintf(stderr,"simulate_ocl: can't open bc file %s\n",bcFileName);
		exit(1);
	}
	fclose(input);
	if(input = fopen(ndFileName,"r")){
		nd = ND_read(input);
	} else {
		fprintf(stderr,"simulate_ocl: can't open nd file %s\n",ndFileName);
		exit(1);
	}
	fclose(input);

	struct ND *res = ND_malloc();
	ND_def_ND(res, nd);
	double *bcv = BCV_malloc(nd->nq);

	double *bcpos = BCPOS_malloc(bc);
	double *bcvel = BCVEL_malloc(bc);
	double *bcrad = BCRAD_malloc(bc);
	long *bcfc = BCFC_malloc(bc);

	set_parameters(bc,nd);
	size_t *ls_item = (size_t *)malloc(2*sizeof(size_t));
	device = create_device_from_file(ls_item, pdFileName); 

	const size_t global_size[2] = {nd->nx,nd->ny};
	const size_t local_size[2] = {ls_item[0],ls_item[1]};

	BCV_def(bc,nd->nq,bcv);
	BCPOS_push(bc,bcpos);
	BCVEL_push(bc,bcvel);
	BCRAD_push(bc,bcrad);
	BCFC_push(bc,bcfc);

//	list_devices();

	fprintf(log_out,"\n[Selected device]: ");
	print_device_name(log_out,device);
	fprintf(log_out,"\n");

	check_workgroup(log_out,global_size,local_size);

	context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
	check_err(err, "Couldn't create a context");
	program = build_program(context, device, programFileName);

	queue = clCreateCommandQueueWithProperties(context, device, 0, &err);
	check_err(err, "Couldn't create a command queue");

	kernel = clCreateKernel(program, KERNEL_FUNC, &err);
	check_err(err, "Couldn't create a kernel");

	print_kernel_info(log_out,kernel, device);
	
	nd_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nd->nq*nd->size*sizeof(double), &nd->m[0], &err);
	res_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, res->nq*res->size*sizeof(double), &res->m[0], &err);
	bcv_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 9*sizeof(double), &bcv[0], &err);
	bcpos_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*bc->nq*sizeof(double), &bcpos[0], &err);
	bcpos_p_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*bc->nq*sizeof(double), &bcpos[0], &err);
	bcvel_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*bc->nq*sizeof(double), &bcvel[0], &err);
	bcrad_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*sizeof(double), &bcrad[0], &err);
	bcfc_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*bc->nq*sizeof(long), &bcfc[0], &err);
	check_err(err, "Couldn't create a buffer");

	err = clSetKernelArg(kernel, 0, sizeof(cl_int), &nd->nq);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_double), &CF);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &nd_buffer);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &res_buffer);
	err |= clSetKernelArg(kernel, 4, sizeof(cl_int), &bc->no);
	err |= clSetKernelArg(kernel, 5, sizeof(cl_int), &bc->nq);
	err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &bcv_buffer);
	err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &bcpos_buffer);
	err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &bcpos_p_buffer);
	err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &bcvel_buffer);
	err |= clSetKernelArg(kernel, 10, sizeof(cl_mem), &bcrad_buffer);
	err |= clSetKernelArg(kernel, 11, sizeof(cl_mem), &bcfc_buffer);
	err |= clSetKernelArg(kernel, 12, sizeof(cl_double), &REFUEL_RTO);
	err |= clSetKernelArg(kernel, 13, sizeof(cl_double), &EAT_RTO);

	check_err(err, "Couldn't create a kernel argument");
	FILE *fp[MP4_NUM];
	char mp4cmd[80];
	double *out_m = (double *)malloc(nd->size*sizeof(double));	
	if(IS_MP4){
		if(check_exist(dirName)){
		for(int i=0;i<MP4_NUM;i++){
			sprintf(mp4cmd,"ffmpeg -y -i - -c:v libx264 -pix_fmt yuv420p %s/%d.mp4 2> /dev/null",dirName,i);
			fp[i] = popen(mp4cmd,"w");
		}
		}else{
			fprintf(stderr,"Output directory not exist\n");
			exit(1);
		}
	}
	err = clEnqueueCopyBuffer(queue,nd_buffer,res_buffer,0,0,nd->nq*nd->size*sizeof(double),0,NULL,NULL);
	check_err(err, "Couldn't copy buffers");
	
	print_parameters(log_out,bc,nd);
	fclose(log_out);
	log_out = fopen(LOG_FILE,"rw");
	if(IS_LOG_PRINT){
		while(getline(&log_line,&log_len,log_out) != EOF){
			printf("%s",log_line);
		}
	}
	if(IS_PROGRESS_PRINT){
		printf("Simulate......");
	}
	fflush(stdout);
	for(int l=0;l<LOOP;l++){
		if(IS_PROGRESS_PRINT){
			printf("\rSimulate......%d/%d",l+1,LOOP);
			fflush(stdout);
		}
		BCFC_init(bcfc,bc->no,bc->nq);
		err = clEnqueueWriteBuffer(queue, bcfc_buffer, CL_TRUE, 0, bc->no*bc->nq*sizeof(long), &bcfc[0], 0, NULL, NULL);
		err = clFinish(queue);
		check_err(err, "Couldn't write the buffer");
		for(int s=0;s<SKP;s++){
			err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, &global_size[0], &local_size[0], 0, NULL,NULL);
			check_err(err, "Couldn't enqueue the kernel");
			err = clFinish(queue);
			err = clEnqueueCopyBuffer(queue,res_buffer,nd_buffer,0,0,nd->nq*nd->size*sizeof(double),0,NULL,NULL);
			check_err(err, "Couldn't copy buffers");
			err = clFinish(queue);
		}
		err = clEnqueueReadBuffer(queue, bcfc_buffer, CL_TRUE, 0, bc->no*bc->nq*sizeof(long), &bcfc[0], 0, NULL, NULL);
		check_err(err, "Couldn't read the buffer");
		err = clFinish(queue);
		
		err = clEnqueueCopyBuffer(queue,bcpos_buffer,bcpos_p_buffer,0,0,bc->no*bc->nq*sizeof(double),0,NULL,NULL);
		check_err(err, "Couldn't copy buffers");
		
		BCFC_pull(bc,bcfc,SKP);
		//ND_probe(debug_out,nd,240,200);
		BC_move_rk4(bc,SKP);
		BC_print(debug_out,bc,0,l*SKP);
		BCPOS_push(bc,bcpos);
		BCVEL_push(bc,bcvel);

		err = clEnqueueWriteBuffer(queue, bcpos_buffer, CL_TRUE, 0, bc->no*bc->nq*sizeof(double), &bcpos[0], 0, NULL, NULL);
		err = clEnqueueWriteBuffer(queue, bcvel_buffer, CL_TRUE, 0, bc->no*bc->nq*sizeof(double), &bcvel[0], 0, NULL, NULL);
		check_err(err, "Couldn't write the buffer");
		err = clFinish(queue);

		if(IS_MP4 || IS_FILE_OUTPUT){
			err = clEnqueueReadBuffer(queue, res_buffer, CL_TRUE, 0, nd->nq*nd->size*sizeof(double), &nd->m[0], 0, NULL, NULL);
			check_err(err, "Couldn't read the buffer");
			err = clFinish(queue);
		}
		if(IS_MP4){
			get_density(out_m,nd);
			nd_ppm_write(out_m,nd->nx,nd->ny,PL_MAX_D,0,1,fp[0]);
			get_ux(out_m,nd);
			nd_ppm_write(out_m,nd->nx,nd->ny,PL_MAX_UX,-1*PL_MAX_UX,1,fp[1]);
			get_uy(out_m,nd);
			nd_ppm_write(out_m,nd->nx,nd->ny,PL_MAX_UY,-1*PL_MAX_UX,1,fp[2]);
		}
		if(IS_FILE_OUTPUT)
		{
		}

	}
	fprintf(log_out,"\rSimulate......completed!! (%dx%d --> %lfs)\n",LOOP,SKP,LOOP*SKP*CT);
	if(IS_PROGRESS_PRINT){
		printf("\rSimulate......completed!! (%dx%d --> %lfs)\n",LOOP,SKP,LOOP*SKP*CT);
	}
	if(IS_SAVE_DATA){
		err = clEnqueueReadBuffer(queue, res_buffer, CL_TRUE, 0, nd->nq*nd->size*sizeof(double), &nd->m[0], 0, NULL, NULL);
		check_err(err, "Couldn't read the buffer");
		err = clFinish(queue);
		sprintf(filename,"%s/fin.nd",dirName);
		if(output = fopen(filename,"w")){
			ND_write(nd,output);
			fclose(output);
		}else{
			fprintf(stderr,"Can't save final ND\n");
		}
	}	

	clReleaseMemObject(nd_buffer);
	clReleaseMemObject(res_buffer);
	clReleaseMemObject(bcpos_buffer);
	clReleaseMemObject(bcvel_buffer);
	clReleaseMemObject(bcrad_buffer);
	clReleaseMemObject(bcfc_buffer);

	clReleaseKernel(kernel);
	clReleaseCommandQueue(queue);
	clReleaseProgram(program);
	clReleaseContext(context);

	fclose(log_out);
	fclose(debug_out);

	ND_free(nd);
	ND_free(res);
	BC_free(bc);
	free(bcpos);
	free(bcvel);
	free(bcrad);
	free(bcfc);
	free(bcv);
}
double msp_get_acc(double dsp,double vel, double ext,double k, double c, double m){
	return (ext - k*dsp - c*vel)/m;
}
void BC_move_rk4(struct BC *bc,double dt){
	struct CY *tmp = NULL;
	double dx,k1,k2,k3,k4;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			dx = tmp->vel[j] * dt;
			k1 = dt*msp_get_acc(tmp->dsp[j],tmp->vel[j],tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
			k2 = dt*msp_get_acc(tmp->dsp[j]+dx/2,tmp->vel[j] + k1/2,tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
			k3 = dt*msp_get_acc(tmp->dsp[j]+dx/2,tmp->vel[j] + k2/2,tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
			k4 = dt*msp_get_acc(tmp->dsp[j]+dx,tmp->vel[j]+k3,tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
			tmp->dsp[j] += dx;
			tmp->pos[j] += dx;
			tmp->vel[j] += (k1+2*k2+2*k3+k4)/6;
			tmp->acc[j] = msp_get_acc(tmp->dsp[j],tmp->vel[j],tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
		}
	}
}
void BC_move_euler(struct BC *bc,double dt){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			tmp->acc[j] = msp_get_acc(tmp->dsp[j],tmp->vel[j],tmp->force[j],tmp->spring,tmp->damp,tmp->mass);
			tmp->vel[j] += tmp->acc[j] * dt;
			tmp->dsp[j] += tmp->vel[j] * dt;
			tmp->pos[j] += tmp->vel[j] * dt;
		}
	}
}
void BC_print(FILE *f,struct BC *bc,int obj,int step){
	struct CY *tmp;
	tmp = (struct CY *)bc->m[obj];
	fprintf(f,"%d %lf ",obj,step*CT);
	fprintf(f,"%lf %lf ",Cforce*tmp->force[0],Cforce*tmp->force[1]);
	fprintf(f,"%lf %lf ",Cacc*tmp->acc[0],Cacc*tmp->acc[1]);
	fprintf(f,"%lf %lf ",Cvel*tmp->vel[0],Cvel*tmp->vel[1]);
	fprintf(f,"%lf %lf ",CL*tmp->dsp[0],CL*tmp->dsp[1]);
	fprintf(f,"%lf %lf\n",CL*tmp->pos[0],CL*tmp->pos[1]);
}
void ND_probe(FILE *f,struct ND *nd, int x, int y){
	double d,ux,uy;
	int idx = x + y*nd->nx;
	for(int vc=0;vc<nd->nq;vc++){
		d += nd->m[vc+idx];
		ux += nd->m[vc+idx]*LC[0+vc*2];
		uy += nd->m[vc+idx]*LC[1+vc*2];
	}
	ux = ux/d;
	uy = uy/d;
	fprintf(f,"%d %d %lf %lf %lf\n",x,y,d,ux,uy);
	fflush(f);
}
double *BCV_malloc(int nq){
	return (double *)malloc(nq*sizeof(double));
}
void BCV_def(struct BC *bc,int nq,double *bcv){
	double v1, v2;
	double cs2 = pow(CS_LTTC,2);
	for (int vc=0;vc<nq;vc++){
		v1 = bc->ux*LC[0+2*vc] + bc->uy*LC[1+2*vc];
		v2 = bc->ux*bc->ux + bc->uy*bc->uy;
		bcv[vc] = W[vc]*bc->dnt*(1 + v1/(cs2) + (v1*v1)/(cs2*cs2*2) - v2/(2*cs2));
	}
	
}
double *BCVEL_malloc(struct BC *bc){
	return (double *)malloc(bc->no*bc->nq*sizeof(double));
}
void BCVEL_push(struct BC *bc, double *bcvel){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			bcvel[j+i*bc->nq] = tmp->vel[j];
		}
	}
}

double *BCPOS_malloc(struct BC *bc){
	return (double *)malloc(bc->no*bc->nq*sizeof(double));
}
void BCPOS_push(struct BC *bc, double *bcpos){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			bcpos[j+i*bc->nq] = tmp->pos[j];
		}
	}
}
long *BCFC_malloc(struct BC *bc){
	return (long *)malloc(bc->no*bc->nq*sizeof(long));
}
void BCFC_push(struct BC *bc, long *bcfc){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			bcfc[j+i*bc->nq] = (long)(tmp->force[j]*FC_OFFSET);
		}
	}
}
void BCFC_init(long *bcfc,int nq,int no){
	for(int i=0;i<no*nq;i++){
			bcfc[i] = 0;
	}
}
void BCFC_pull(struct BC *bc,long *bcfc, double dt){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			tmp->force[j] = (double)(bcfc[j+i*bc->nq])/(dt*FC_OFFSET);
		}
	}
}
double *BCRAD_malloc(struct BC *bc){
	return (double *)malloc(bc->no*sizeof(double));
}
double *BCRAD_push(struct BC *bc,double *bcrad){
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		bcrad[i] = tmp->rad;
	}
}

void print_device_name(FILE *f, cl_device_id device){
	char* value;
	size_t valueSize;
	// print device name
	clGetDeviceInfo(device, CL_DEVICE_NAME, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_NAME, valueSize, value, NULL);
	fprintf(f,"%s", value);
	free(value);
}
void print_kernel_info(FILE *f,cl_kernel kernel, cl_device_id device){
	size_t kernel_value, err;
	fprintf(f,"[Kernel info]:\n");
	err = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), NULL, &kernel_value);
	check_err(err, "Couldn't check kernel value");
	fprintf(f,"\tCL_KERNEL_WORK_GROUP_SIZE: %d\n",kernel_value);
}

void list_devices(){

	char* value;
	size_t valueSize;
	cl_int platform_count;
	cl_platform_id* platforms;
	cl_int device_count;
	cl_device_id *devices;
	cl_device_id res_dev;
	cl_int max_compute_units;
	int err;
	size_t *val;

	/* Identify all platform */
	err = clGetPlatformIDs(0, NULL, &platform_count);
	platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platform_count);
	clGetPlatformIDs(platform_count, platforms, NULL);
	check_err(err,"Couldn't identigy a platform");

	printf("[Opencl devices info]:\n");	
	for (int i=0;i<platform_count;i++){
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &device_count);
		devices = (cl_device_id*) malloc(sizeof(cl_device_id) * device_count);
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, device_count, devices, NULL);
		for(int j=0;j<device_count;j++){
			printf("\t");
			print_device_name(stdout,devices[j]);
			printf("\n");
			// print device CL_DEVICE_MAX_WORK_ITEM_SIZES
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_ITEM_SIZES, 0, NULL, &valueSize);
			val = (size_t*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_ITEM_SIZES, valueSize, val, NULL);
			printf("\t\tCL_DEVICE_MAX_WORK_ITEM_SIZES:  %d/%d/%d\n", val[0], val[1], val[2]);
			free(val);
			// print device CL_DEVICE_MAX_WORK_GROUP_SIZE
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, 0, NULL, &valueSize);
			val = (size_t*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, valueSize, &val[0], NULL);
			printf("\t\tCL_DEVICE_MAX_WORK_GROUP_SIZE:  %d\n", val[0]);
			free(val);
			// print device CL_DEVICE_MAX_COMPUTE_UNITS
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, 0, NULL, &valueSize);
			val = (size_t*) malloc(valueSize);
			clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, valueSize, &val[0], NULL);
			printf("\t\tCL_DEVICE_MAX_COMPUTE_UNITS:  %d\n", val[0]);
		}
	}
}

void check_workgroup(FILE *f,const size_t *gs,const size_t *ls){
	fprintf(f,"[Workgroup info]\n");
	fprintf(f,"\tglobal_size: %d/%d (total: %d)\n",gs[0],gs[1],gs[0]*gs[1]);
	fprintf(f,"\tlocal_size: %d/%d (total: %d)\n",ls[0],ls[1],ls[0]*ls[1]);
	if(gs[0] % ls[0] != 0){
		fprintf(stderr,"D 1 GS not divided perfectly by LS\n");
		exit(1);
	}
	if(gs[1] % ls[1] != 0){
		fprintf(stderr,"D 2 GS not divided perfectly by LS\n");
		exit(1);
	}
	fprintf(f,"\tworkgroups: %d/%d (total: %d)\n",gs[0]/ls[0],gs[1]/ls[1],gs[0]*gs[1]/ls[0]/ls[1]);
}

cl_device_id create_device_from_file(size_t *ls, char* file){
	FILE *f;
	size_t p, d, l1, l2;
	cl_device_id device;
	f = fopen(file, "r");
	if (f) {
		fscanf(f,"%d %d %d %d",&p,&d,&l1,&l2);
		ls[0] = l1;
		ls[1] = l2;
		device = create_device(p,d);
	}else{
		printf("can't open pdFile %s\n",file);
		exit(1);
	}
	fclose(f);
	return device;
}

cl_device_id create_device(int sel_plat, int sel_dev) {

	char* value;
	size_t valueSize;
	cl_int platform_count;
	cl_platform_id* platforms;
	cl_int device_count;
	cl_device_id *devices;
	cl_device_id res_dev;
	cl_int max_compute_units;
	int err;
	size_t *val;

	/* Identify all platform */
	err = clGetPlatformIDs(0, NULL, &platform_count);
	platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platform_count);
	clGetPlatformIDs(platform_count, platforms, NULL);
	check_err(err,"Couldn't identigy a platform");

	if(sel_plat+1 > platform_count){
		fprintf(stderr, "selected platform not exist: p%d\n",sel_plat+1);
		exit(1);
	}

	for (int i=0;i<platform_count;i++){
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &device_count);
		devices = (cl_device_id*) malloc(sizeof(cl_device_id) * device_count);
		clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, device_count, devices, NULL);
		for(int j=0;j<device_count;j++){
			if(i==sel_plat){
				if(sel_dev+1 > device_count){
					fprintf(stderr, "selected device not exist: d%d\n",sel_dev+1);
					exit(1);
				}
				if(j==sel_dev){
					res_dev = devices[j];
				}
			}
		}
	}
	if (res_dev == NULL){
		fprintf(stderr, "Selected devices (P%d D%d) not found.",sel_plat,sel_dev);
	}
	return res_dev;
}

/* Create program from a file and compile it */
cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename) {

	cl_program program;
	FILE *program_handle;
	char *program_buffer, *program_log;
	size_t program_size, log_size;
	int err;

	/* Read program file and place content into buffer */
	program_handle = fopen(filename, "r");
	if(program_handle == NULL) {
		perror("Couldn't find the program file");
		exit(1);
	}
	fseek(program_handle, 0, SEEK_END);
	program_size = ftell(program_handle);
	rewind(program_handle);
	program_buffer = (char*)malloc(program_size + 1);
	program_buffer[program_size] = '\0';
	fread(program_buffer, sizeof(char), program_size, program_handle);
	fclose(program_handle);

	/* Create program from file 

	   Creates a program from the source code in the add_numbers.cl file. 
	   Specifically, the code reads the file's content into a char array 
	   called program_buffer, and then calls clCreateProgramWithSource.
	   */
	program = clCreateProgramWithSource(ctx, 1, 
			(const char**)&program_buffer, &program_size, &err);
	if(err < 0) {
		perror("Couldn't create the program");
		exit(1);
	}
	free(program_buffer);

	/* Build program 

	   The fourth parameter accepts options that configure the compilation. 
	   These are similar to the flags used by gcc. For example, you can 
	   define a macro with the option -DMACRO=VALUE and turn off optimization 
	   with -cl-opt-disable.
	   */
	err = clBuildProgram(program, 0, NULL, "-I. -I /usr/include -I /usr/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include ", NULL, NULL);
	if(err < 0) {

		/* Find size of log and print to std output */
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		program_log = (char*) malloc(log_size + 1);
		program_log[log_size] = '\0';
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 
				log_size + 1, program_log, NULL);
		printf("%s\n", program_log);
		free(program_log);
		exit(1);
	}

	return program;
}
const char *getErrorString(cl_int error)
{
	switch(error){
		// run-time and JIT compiler errors
		case 0: return "CL_SUCCESS";
		case -1: return "CL_DEVICE_NOT_FOUND";
		case -2: return "CL_DEVICE_NOT_AVAILABLE";
		case -3: return "CL_COMPILER_NOT_AVAILABLE";
		case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case -5: return "CL_OUT_OF_RESOURCES";
		case -6: return "CL_OUT_OF_HOST_MEMORY";
		case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case -8: return "CL_MEM_COPY_OVERLAP";
		case -9: return "CL_IMAGE_FORMAT_MISMATCH";
		case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case -11: return "CL_BUILD_PROGRAM_FAILURE";
		case -12: return "CL_MAP_FAILURE";
		case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case -15: return "CL_COMPILE_PROGRAM_FAILURE";
		case -16: return "CL_LINKER_NOT_AVAILABLE";
		case -17: return "CL_LINK_PROGRAM_FAILURE";
		case -18: return "CL_DEVICE_PARTITION_FAILED";
		case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

			  // compile-time errors
		case -30: return "CL_INVALID_VALUE";
		case -31: return "CL_INVALID_DEVICE_TYPE";
		case -32: return "CL_INVALID_PLATFORM";
		case -33: return "CL_INVALID_DEVICE";
		case -34: return "CL_INVALID_CONTEXT";
		case -35: return "CL_INVALID_QUEUE_PROPERTIES";
		case -36: return "CL_INVALID_COMMAND_QUEUE";
		case -37: return "CL_INVALID_HOST_PTR";
		case -38: return "CL_INVALID_MEM_OBJECT";
		case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case -40: return "CL_INVALID_IMAGE_SIZE";
		case -41: return "CL_INVALID_SAMPLER";
		case -42: return "CL_INVALID_BINARY";
		case -43: return "CL_INVALID_BUILD_OPTIONS";
		case -44: return "CL_INVALID_PROGRAM";
		case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
		case -46: return "CL_INVALID_KERNEL_NAME";
		case -47: return "CL_INVALID_KERNEL_DEFINITION";
		case -48: return "CL_INVALID_KERNEL";
		case -49: return "CL_INVALID_ARG_INDEX";
		case -50: return "CL_INVALID_ARG_VALUE";
		case -51: return "CL_INVALID_ARG_SIZE";
		case -52: return "CL_INVALID_KERNEL_ARGS";
		case -53: return "CL_INVALID_WORK_DIMENSION";
		case -54: return "CL_INVALID_WORK_GROUP_SIZE";
		case -55: return "CL_INVALID_WORK_ITEM_SIZE";
		case -56: return "CL_INVALID_GLOBAL_OFFSET";
		case -57: return "CL_INVALID_EVENT_WAIT_LIST";
		case -58: return "CL_INVALID_EVENT";
		case -59: return "CL_INVALID_OPERATION";
		case -60: return "CL_INVALID_GL_OBJECT";
		case -61: return "CL_INVALID_BUFFER_SIZE";
		case -62: return "CL_INVALID_MIP_LEVEL";
		case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
		case -64: return "CL_INVALID_PROPERTY";
		case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
		case -66: return "CL_INVALID_COMPILER_OPTIONS";
		case -67: return "CL_INVALID_LINKER_OPTIONS";
		case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

			  // extension errors
		case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
		case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
		case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
		case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
		case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
		case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
		default: return "Unknown OpenCL error";
	}
}
void check_err(int err, char *str){
	if (err < 0 ){
		fprintf(stderr,"%s: %s\n",str,getErrorString(err));
		exit(1);
	}	
}
