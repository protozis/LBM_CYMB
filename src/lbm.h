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
struct GP {
	int no;
	void **obj;
}
struct CY {
	int nq;
	double rist;
	double damp;
	double mass;
	int32_t *force;
	double *acc;
	double *vel;
	double *dsp;
	double *pos;
};
struct GP *GP_malloc();

void CY_def(struct CY *cy,int no,int nq):
void CY_free(struct CY *cy);
void CY_init(struct CY *cy, double force, double acc, double vel, double dsp);
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
