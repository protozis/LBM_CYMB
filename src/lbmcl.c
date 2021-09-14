#define KERNEL_FUNC "propagate"
#define CL_TARGET_OPENCL_VERSION 300
#define CY_PAR_NUM 4
#define CY_KIE_NUM 5

#ifdef MAC
#include<OpenCL/cl.h>
#else
#include<CL/cl.h>
#endif

#include"lbm.h"
#include"lbmcl.h"
#include<stdio.h>
#include<string.h>

//Default
size_t LOOP = 1;
size_t SKP = 1;
double CF = 1;
double SL = 1;
int IS_MP4 = 0;
int IS_SAVE_DATA = 0;
int IS_FILE_OUTPUT = 0;
char ND_FILE[80];
char BC_FILE[80];
char PD_FILE[80];
char DIR_NAME[80];
char PROGRAM_FILE[80];

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
			}else if(strcmp(par,"SL") == 0){
				SL = strtod(val,&ptr);
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
			}else if(strcmp(par,"DIR_NAME") == 0){
				 memcpy(DIR_NAME,val,strlen(val)+1);
			}else if(strcmp(par,"PROGRAM_FILE") == 0){
				 memcpy(PROGRAM_FILE,val,strlen(val)+1);
			}
		}

	}else{
		fprintf(stderr,"read_config: can't open conf file\n");
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
	cl_mem bcp_buffer;
	cl_mem bck_buffer;

	struct BC *bc = NULL;
	struct ND *nd = NULL;
	FILE *output;
	FILE *input;
	char filename[80];
	if(input = fopen(bcFileName,"r")){
		bc = BC_read(input);
		fclose(input);
	} else {
		fprintf(stderr,"simulate_ocl: can't open bc file\n");
		exit(0);
	}
	if(input = fopen(ndFileName,"r")){
		nd = ND_read(input);
	} else {
		fprintf(stderr,"simulate_ocl: can't open bc file\n");
		exit(0);
	}


	struct ND *res = ND_malloc();
	ND_def_ND(res, nd);
	double *bcv = BCV_malloc(nd->nq);
	double *bcp = BCP_malloc(bc);
	double *bck = BCK_malloc(bc);
	BCKP_def(bc,bcp,bck);
	BCV_def(bc,nd->nq,bcv);

	size_t *ls_item = (size_t *)malloc(2*sizeof(size_t));
	device = create_device_from_file(ls_item, pdFileName); 

	const size_t global_size[2] = {nd->nx,nd->ny};
	const size_t local_size[2] = {ls_item[0],ls_item[1]};
	
//	list_devices();

	printf("\n[Selected device]: ");
	print_device_name(device);
	printf("\n");

	check_workgroup(global_size,local_size);

	context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
	check_err(err, "Couldn't create a context");
	program = build_program(context, device, programFileName);

	queue = clCreateCommandQueueWithProperties(context, device, 0, &err);
	check_err(err, "Couldn't create a command queue");

	kernel = clCreateKernel(program, KERNEL_FUNC, &err);
	check_err(err, "Couldn't create a kernel");

	print_kernel_info(kernel, device);

	nd_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nd->nq*nd->size*sizeof(double), &nd->m[0], &err);
	res_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, res->nq*res->size*sizeof(double), &res->m[0], &err);
	bcv_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 9*sizeof(double), &bcv[0], &err);
	bcp_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*CY_PAR_NUM*sizeof(double), &bcp[0], &err);
	bck_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, bc->no*bc->nq*CY_KIE_NUM*sizeof(double), &bck[0], &err);
	check_err(err, "Couldn't create a buffer");

	err = clSetKernelArg(kernel, 0, sizeof(cl_uint), &nd->nq);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_double), &SL);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_double), &CF);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &nd_buffer);
	err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &res_buffer);
	err |= clSetKernelArg(kernel, 5, sizeof(cl_uint), &bc->no);
	err |= clSetKernelArg(kernel, 6, sizeof(cl_uint), &bc->nq);
	err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &bcv_buffer);
	err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &bcp_buffer);
	err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &bck_buffer);

	check_err(err, "Couldn;t create a kernel argument");
	char mp4cmd[80];
	if(IS_MP4){
	}

	printf("[Parameters]:\n");

	printf("Simulate......");
	fflush(stdout);
	for(int l=0;l<LOOP;l++){
		check_err(err, "Couldn't write the buffer");

		printf("\rSimulate......%d/%d",l+1,LOOP);
		fflush(stdout);

		for(int s=0;s<SKP;s++){

			err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, &global_size[0], &local_size[0], 0, NULL,NULL);
			check_err(err, "Couldn't enqueue the kernel");

			err = clFinish(queue);
			check_err(err, "Queue not finish");

		}
		err = clEnqueueReadBuffer(queue, res_buffer, CL_TRUE, 0, nd->nq*nd->size*sizeof(double), &nd->m[0], 0, NULL, NULL);
		err = clEnqueueReadBuffer(queue, bck_buffer, CL_TRUE, 0, bc->no*bc->nq*CY_KIE_NUM*sizeof(double), &bck[0], 0, NULL, NULL);
		check_err(err, "Couldn't read the buffer");

		err = clFinish(queue);
		check_err(err, "Queue not finish");
		if(IS_MP4){
		}
		if(IS_FILE_OUTPUT)
		{
			sprintf(filename,"%s/%.4d",dirName,l);
			output = fopen(filename,"w");
			ND_write(nd,output);
			fclose(output);
		}

	}
	printf("\rSimulate......completed!! (x%d)\n",LOOP);

	if(IS_SAVE_DATA){
		sprintf(filename,"%s/fin.nd",dirName);
		output = fopen(filename,"w");
		ND_write(nd,output);
		fclose(output);
	}	

	clReleaseMemObject(nd_buffer);
	clReleaseMemObject(res_buffer);
	clReleaseMemObject(bcp_buffer);
	clReleaseMemObject(bck_buffer);

	clReleaseKernel(kernel);
	clReleaseCommandQueue(queue);
	clReleaseProgram(program);
	clReleaseContext(context);

	ND_free(nd);
	ND_free(res);
	BC_free(bc);
	free(bcp);
	free(bck);
}
double *BCV_malloc(uint nq){
	return (double *)malloc(nq*sizeof(double));
}
double *BCK_malloc(struct BC *bc){
	return (double *)malloc(bc->no*bc->nq*CY_KIE_NUM*sizeof(double));
}
double *BCP_malloc(struct BC *bc){
	return (double *)malloc(bc->no*CY_PAR_NUM*sizeof(double));
}
void BCV_def(struct BC *bc,uint nq,double *bcv){
	double v1, v2;
	double cs2 = SL*SL/3;
	for (int vc=0;vc<nq;vc++){
		v1 = bc->ux*LC[0+2*vc] + bc->uy*LC[1+2*vc];
		v2 = bc->ux*bc->ux + bc->uy*bc->uy;
		bcv[vc] = W[vc]*bc->dnt*(1 + v1/(cs2) + (v1*v1)/(cs2*cs2*2) - v2/(2*cs2));
	}
	
}
void BCKP_def(struct BC *bc, double *bcp, double *bck){
	int pn = CY_PAR_NUM;
	int kn = CY_KIE_NUM*bc->nq;
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		bcp[pn*i+0] = tmp->rist;
		bcp[pn*i+1] = tmp->damp;
		bcp[pn*i+2] = tmp->mass;
		bcp[pn*i+3] = tmp->rad;
		for(int j=0;j<bc->nq;j++){
			bck[kn*i+0*bc->nq+j] = tmp->force[j];
			bck[kn*i+1*bc->nq+j] = tmp->acc[j];
			bck[kn*i+2*bc->nq+j] = tmp->vel[j];
			bck[kn*i+3*bc->nq+j] = tmp->dsp[j];
			bck[kn*i+4*bc->nq+j] = tmp->pos[j];
		}
	}
}
void BCKP_update(struct BC *bc, double *bcp, double *bck){
	int pn = CY_PAR_NUM;
	int kn = CY_KIE_NUM*bc->nq;
	struct CY *tmp = NULL;
	for(int i=0;i<bc->no;i++){
		tmp = ((struct CY *)bc->m[i]);
		for(int j=0;j<bc->nq;j++){
			tmp->force[j] = bck[kn*i+0*bc->nq+j];
			tmp->acc[j] = bck[kn*i+1*bc->nq+j];
			tmp->vel[j] = bck[kn*i+2*bc->nq+j];
			tmp->dsp[j] = bck[kn*i+3*bc->nq+j];
			tmp->pos[j] = bck[kn*i+4*bc->nq+j];
		}
	}
}

void print_device_name( cl_device_id device){
	char* value;
	size_t valueSize;
	// print device name
	clGetDeviceInfo(device, CL_DEVICE_NAME, 0, NULL, &valueSize);
	value = (char*) malloc(valueSize);
	clGetDeviceInfo(device, CL_DEVICE_NAME, valueSize, value, NULL);
	printf("%s", value);
	free(value);
}
void print_kernel_info(cl_kernel kernel, cl_device_id device){
	size_t kernel_value, err;
	printf("[Kernel info]:\n");
	err = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), NULL, &kernel_value);
	check_err(err, "Couldn't check kernel value");
	printf("\tCL_KERNEL_WORK_GROUP_SIZE: %d\n",kernel_value);
}

void list_devices(){

	char* value;
	size_t valueSize;
	cl_uint platform_count;
	cl_platform_id* platforms;
	cl_uint device_count;
	cl_device_id *devices;
	cl_device_id res_dev;
	cl_uint max_compute_units;
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
			print_device_name(devices[j]);
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

void check_workgroup(const size_t *gs,const size_t *ls){
	printf("[Workgroup info]\n");
	printf("\tglobal_size: %d/%d (total: %d)\n",gs[0],gs[1],gs[0]*gs[1]);
	printf("\tlocal_size: %d/%d (total: %d)\n",ls[0],ls[1],ls[0]*ls[1]);
	if(gs[0] % ls[0] != 0){
		fprintf(stderr,"D 1 GS not divided perfectly by LS\n");
		exit(1);
	}
	if(gs[1] % ls[1] != 0){
		fprintf(stderr,"D 2 GS not divided perfectly by LS\n");
		exit(1);
	}
	printf("\tworkgroups: %d/%d (total: %d)\n",gs[0]/ls[0],gs[1]/ls[1],gs[0]*gs[1]/ls[0]/ls[1]);
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
		printf("can't open pdFile\n");
		exit(1);
	}
	fclose(f);
	return device;
}

cl_device_id create_device(int sel_plat, int sel_dev) {

	char* value;
	size_t valueSize;
	cl_uint platform_count;
	cl_platform_id* platforms;
	cl_uint device_count;
	cl_device_id *devices;
	cl_device_id res_dev;
	cl_uint max_compute_units;
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
