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

struct CY {
	double *pos;
	double rist;
	double damp;
	double mass;
	int32_t *force;
	double *acc;
	double *vel;
	double *dsp;
};

void ND_propagate(struct ND *output, struct ND *input, struct BC *bc, struct FC *fc, double sl, double cf);
void ND_eq(struct ND *eq , struct ND *nd, double sl);
void ND_collision(struct ND *res, struct ND *old, struct ND *eq, struct BC *bc, double cf);
void ND_streaming(struct ND *res, struct ND *old, struct BC *bc, struct FC *fc);
void ND_free(struct ND *nd);
void ND_init(double *tmp, double D, double ux, double uy, double sl);
void ND_def_ND(struct ND *nd, struct ND *ndp);
void ND_copy(struct ND *nd,struct ND *ref);
void ND_def(struct ND *nd, int nx, int ny, int nq);
struct ND *ND_malloc();

void get_density(double *tmp, struct ND * nd);
void get_ux(double *tmp, struct ND * nd);
void get_uy(double *tmp, struct ND * nd);

int isnum(char c);

void jet_colormap(double num, double max, int *c);


struct ND *ND_read(FILE *f);
void ND_write(struct ND *nd, FILE *f);


void nd_pgm_write(double *m, int nx, int ny, int pmax, int scale, FILE *f);
void nd_ppm_write(double *m, int nx, int ny, double pmax, int scale, FILE *f);

#endif
