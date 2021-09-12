#include<stdlib.h>
#include<stdint.h>
#include<string.h>
#include<regex.h>
#include "lbm.h"

const double LC[18] = {
	-1,-1,
	0,-1,
	1,-1,
	-1,0,
	0,0,
	1,0,
	-1,1,
	0,1,
	1,1
};

const double W[9] ={ 
	0.027777,
	0.111111,
	0.027777,
	0.111111,
	0.444444,
	0.111111,
	0.027777,
	0.111111,
	0.027777
};
struct BC *BC_read(FILE *f){
	struct BC *bc = BC_malloc();
	char flag[80];
	uint cnt = 0;
	int c;
	fscanf(f,"bc_no %d\n",&bc->no);
	fscanf(f,"bc_nq %d\n",&bc->nq);
	BC_def(bc,bc->no);
	while(fscanf(f,"%s {\n",flag) >= 0){
		if(strcmp(flag,"BCV") == 0){
			fscanf(f,"dnt %lf\n",&bc->dnt);
			fscanf(f,"ux %lf\n",&bc->ux);
			fscanf(f,"uy %lf\n",&bc->uy);
		} else if (strcmp(flag,"CY") == 0){
			bc->m[cnt] = CY_read(f);
			cnt++;
		}
		fscanf(f,"}\n");
	}
	return bc;
}
struct BC *BC_malloc(){
	return (struct BC *)malloc(sizeof(struct BC));
}
void BC_def(struct BC *bc,int no){
	bc->m = (void **)malloc(no*sizeof(void *));
}
void BC_free(struct BC *bc){
	free(bc->m);
	free(bc);
}
struct CY *CY_read(FILE *f){
	struct CY *cy = CY_malloc();
	CY_def(cy);
	fscanf(f,"rist %lf\n",&cy->rist);
	fscanf(f,"damp %lf\n",&cy->damp);
	fscanf(f,"mass %lf\n",&cy->mass);
	fscanf(f,"rad %lf\n",&cy->rad);
	fscanf(f,"force %lf %lf\n",&cy->force[0],&cy->force[1]);
	fscanf(f,"acc %lf %lf\n",&cy->acc[0],&cy->acc[1]);
	fscanf(f,"vel %lf %lf\n",&cy->vel[0],&cy->vel[1]);
	fscanf(f,"dsp %lf %lf\n",&cy->dsp[0],&cy->dsp[1]);
	fscanf(f,"pos %lf %lf\n",&cy->pos[0],&cy->pos[1]);
	return cy;
}
void CY_init(struct CY *cy, double force, double acc, double vel, double dsp){
	for(int i=0;i<2;i++){
		cy->force[i] = force;
		cy->acc[i] = acc;
		cy->vel[i] = vel;
		cy->dsp[i] = dsp;
	}
}
void CY_def(struct CY *cy){
	cy->force = (double *)malloc(2*sizeof(double));
	cy->acc = (double *)malloc(2*sizeof(double));
	cy->vel = (double *)malloc(2*sizeof(double));
	cy->dsp = (double *)malloc(2*sizeof(double));
	cy->pos = (double *)malloc(2*sizeof(double));
}
void CY_free(struct CY *cy){
	free(cy->pos);
	free(cy->force);
	free(cy->acc);
	free(cy->vel);
	free(cy->dsp);
	free(cy);
}
struct CY *CY_malloc(){
	return (struct CY *)malloc(sizeof(struct CY));
}

void get_density(double *tmp, struct ND * nd){
	for(int i=0;i<nd->size;i++){
		tmp[i] = 0;
		for(int j=0;j<nd->nq;j++){
			tmp[i] = tmp[i] + nd->m[j+i*nd->nq];
		}
	}
}

void get_ux(double *tmp, struct ND * nd){
	double *D = (double*)malloc(nd->size*sizeof(double));
	get_density(D,nd); 
	for(int i=0;i<nd->size;i++){
		tmp[i] = 0;
		for(int j=0;j<nd->nq;j++){
			tmp[i] = tmp[i] + nd->m[j+i*nd->nq]*LC[0+2*j];
		}
		if(D[i] != 0){
			tmp[i] = tmp[i]/D[i];
		} else {
			tmp[i] = 0;
		}
	}
	free(D);
}

void get_uy(double *tmp, struct ND * nd){
	double *D = (double*)malloc(nd->size*sizeof(double));
	get_density(D,nd);
	for(int i=0;i<nd->size;i++){
		tmp[i] = 0;
		for(int j=0;j<nd->nq;j++){
			tmp[i] = tmp[i] + nd->m[j+i*nd->nq]*LC[1+2*j];
		}
		if(D[i] != 0){
			tmp[i] = tmp[i]/D[i];
		} else {
			tmp[i] = 0;
		}
	}
	free(D);
}

void ND_init(double *tmp, double D, double ux, double uy, double sl){
	double d1,d2;
	double cs2 = sl*sl/3;
	d2 = ux*ux + uy*uy;
	for(int i=0;i<9;i++){
		d1 = ux*LC[0+2*i]+uy*LC[1+2*i];
		tmp[i]=W[i]*D*(1+d1/(cs2)+(d1*d1)/(cs2*cs2*2)-d2/(2*cs2));
	}
}

void ND_free(struct ND *nd){
	free(nd->m);
	free(nd);
}

int isnum(char c){
	int tmp = 0;
	if( ((int)c-48) >= 0 && ((int)c-48) <= 9 ){
		tmp = 1;
	}
	return tmp;
}

struct ND  *ND_malloc(){
	return (struct ND *)malloc(sizeof(struct ND));
}
void ND_def_ND(struct ND *nd, struct ND *ndp){
	nd->nx = ndp->nx;
	nd->ny = ndp->ny;
	nd->nq = ndp->nq;
	nd->size = ndp->size;
	nd->m = (double *)malloc(nd->size*nd->nq*sizeof(double));
}
void ND_copy(struct ND *nd,struct ND *ref){
	for(int i=0;i<ref->nq*ref->size;i++){
		nd->m[i] = ref->m[i];
	}
}
void ND_def(struct ND *nd, int nx, int ny, int nq){
	nd->nx = nx;
	nd->ny = ny;
	nd->nq = nq;
	nd->size = nx*ny;
	nd->m = (double *)malloc(nd->size*nd->nq*sizeof(double));
}
struct ND *ND_read(FILE *f){
	int nx, ny, nq;
	fscanf(f, "%d",&nx);
	fscanf(f, "%d",&ny);
	fscanf(f, "%d",&nq);
	struct ND *ndf = ND_malloc();
	ND_def(ndf,nx,ny,nq);
	for(int i=0;i < (ndf->size*ndf->nq);i++){
		fscanf(f,"%lf",&ndf->m[i]);
	}
	fclose(f);
	return ndf;
}

void ND_write(struct ND *nd, FILE *f){
	fprintf(f,"%d %d %d\n",nd->nx,nd->ny,nd->nq);
	for(int i=0;i< nd->size;i++){
		for(int j=0;j<nd->nq;j++){
			fprintf(f,"%lf ",nd->m[j+i*nd->nq]);
		}
		fprintf(f,"\n");
	}
}
