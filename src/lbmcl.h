#ifndef LBMCL
#define LBMCL

#define CL_TARGET_OPENCL_VERSION 300

#ifdef MAC
#include<OpenCL/cl.h>
#else
#include<CL/cl.h>
#endif

#include<stdio.h>
extern size_t LOOP;
extern size_t SKP;
extern double CF;
extern double SL;
extern int IS_MP4;
extern int IS_SAVE_DATA;
extern int IS_FILE_OUTPUT;
extern char ND_FILE[80];
extern char BC_FILE[80];
extern char PD_FILE[80];
extern char DIR_NAME[80];
extern char PROGRAM_FILE[80];

void update_config(char* filename);
void simulate_ocl(char* ndFileName, char* bcFileName, char* pdFileName, char* dirName, char* programFileName);
double *BCV_malloc();
double *BCK_malloc(struct BC *bc);
double *BCP_malloc(struct BC *bc);
void BCV_def(struct BC *bc,double *bcv);
void BCKP_def(struct BC *bc, double *bcp, double *bck);
void BCKP_update(struct BC *bc, double *bcp, double *bck);
void list_devices();
void print_device_name(cl_device_id device);
void print_kernel_info(cl_kernel kernel, cl_device_id device);
void check_workgroup(const size_t *gs, const size_t *ls);
cl_device_id create_device_from_file(size_t * ls, char* file);
cl_device_id create_device(int sel_plat, int sel_dev);

cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename);

const char *getErrorString(cl_int error);

void check_err(int err, char *str);
#endif
