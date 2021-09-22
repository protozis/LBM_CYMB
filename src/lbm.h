#ifndef LBM
#define LBM

#include<stdio.h>

extern const double LC[18];

extern const double W[9];

struct ND {
	int nx;
	int ny;
	int nq;
	int size;
	double *m;
};
struct BC {
	int no;
	int nq;
	double dnt;
	double ux;
	double uy;
	void **m;
};
struct CY {
	double spring;
	double damp;
	double mass;
	double rad;
	double *force;
	double *acc;
	double *vel;
	double *dsp;
	double *pos;
};
struct BC *BC_read(FILE *f);
struct BC *BC_malloc();
void BC_def(struct BC *bc,int no);
void BC_free(struct BC *bc);

void CY_def(struct CY *cy);
void CY_free(struct CY *cy);
void CY_init(struct CY *cy, double force, double acc, double vel, double dsp);
struct CY *CY_read(FILE *f);
struct CY *CY_malloc();

void ND_free(struct ND *nd);
void ND_init(double *tmp, double D, double ux, double uy, double sl);
void ND_def_ND(struct ND *nd, struct ND *ndp);
void ND_copy(struct ND *nd,struct ND *ref);
void ND_def(struct ND *nd, int nx, int ny, int nq);
struct ND *ND_malloc();
struct ND *ND_read(FILE *f);
void ND_write(struct ND *nd, FILE *f);

void get_density(double *tmp, struct ND * nd);
void get_ux(double *tmp, struct ND * nd);
void get_uy(double *tmp, struct ND * nd);

int isnum(char c);

#endif
