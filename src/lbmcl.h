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
extern double CS;
extern double CL;
extern double CT;
extern double CD;
extern double MA;
extern int IS_MP4;
extern int IS_SAVE_DATA;
extern int IS_FILE_OUTPUT;
extern char ND_FILE[80];
extern char BC_FILE[80];
extern char PD_FILE[80];
extern char OUTPUT_DIR[80];
extern char PROGRAM_FILE[80];
extern double REFUEL_RTO;
extern double EAT_RTO;
extern char LOG_FILE[80];
extern char DEBUG_FILE[80];
extern int IS_LOG_PRINT;
extern int IS_PROGRESS_PRINT;

void update_config(char* filename);
void set_parameters(struct BC *bc, struct ND *nd);
void print_parameters(FILE *f,struct BC *bc,struct ND *nd);
void simulate_ocl(char* ndFileName, char* bcFileName, char* pdFileName, char* dirName, char* programFileName);
void ND_probe(FILE *f,struct ND *nd, int x, int y);
double msp_get_acc(double dsp,double vel, double ext,double k, double c, double m);
void BC_move_rk4(struct BC *bc,double dt);
void BC_move_euler(struct BC *bc,double dt);
void BC_print(FILE *f,struct BC *bc,int obj,int step);
double *BCV_malloc(int nq);
void BCV_def(struct BC *bc, int nq, double *bcv);
double *BCVEL_malloc(struct BC *bc);
void BCVEL_push(struct BC *bc, double *bcvel);
double *BCPOS_malloc(struct BC *bc);
void BCPOS_push(struct BC *bc, double *bcpos);
long *BCFC_malloc(struct BC *bc);
void BCFC_push(struct BC *bc, long *bcfc);
void BCFC_init(long *bcfc,int nq,int no);
void BCFC_pull(struct BC *bc,long *bcfc,double dt);
double *BCRAD_malloc(struct BC *bc);
double *BCRAD_push(struct BC *bc,double *bcrad);
void nd_ppm_write(double *m, int nx, int ny, double pmax, int scale, FILE *f);
void jet_colormap(double num, double min, double max,int *c);
void list_devices();
void print_device_name(FILE *f,cl_device_id device);
void print_kernel_info(FILE *f,cl_kernel kernel, cl_device_id device);
void check_workgroup(FILE *f,const size_t *gs, const size_t *ls);
cl_device_id create_device_from_file(size_t * ls, char* file);
cl_device_id create_device(int sel_plat, int sel_dev);

cl_program build_program(cl_context ctx, cl_device_id dev, const char* filename);

const char *getErrorString(cl_int error);

void check_err(int err, char *str);
#endif
