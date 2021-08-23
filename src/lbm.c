#include<stdlib.h>
#include<stdint.h>
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

void CY_def(struct CY *cy,int dims){
	cy->pos = (double *)malloc(dims*sizeof(double));
	cy->force = (int32_t *)malloc(dims*sizeof(int32_t));
	cy->acc = (double *)malloc(dims*sizeof(double));
	cy->vel = (double *)malloc(dims*sizeof(double));
	cy->dsp = (double *)malloc(dims*sizeof(double));
}
void CY_free(struct CY *cy){
	free(cy->pos);
	free(cy->force);
	free(cy->acc);
	free(cy->vel);
	free(cy->dep);
	free(cy);
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

void ND_propagate(struct ND *output, struct ND *input, struct BC *bc, struct FC *fc, double sl, double cf){
	struct ND *ns = ND_malloc();
	struct ND *eq = ND_malloc();
	ND_def_ND(ns,input);
	ND_def_ND(eq,input);
	ND_streaming(ns, input, bc, fc);
	ND_eq(eq,ns,sl);
	ND_collision(output, ns, eq, bc, cf);
	ND_free(ns);
	ND_free(eq);
}

void ND_eq(struct ND *eq, struct ND *nd, double sl){
	double *D;
	double *D2;
	double *ux;
	double *uy;
	D = (double *)malloc(nd->size*sizeof(double));
	D2 = (double *)malloc(nd->size*sizeof(double));
	ux = (double *)malloc(nd->size*sizeof(double));
	uy = (double *)malloc(nd->size*sizeof(double));
	get_density(D, nd);
	get_ux(ux, nd);
	get_uy(uy, nd);
	int idx;
	double d1,d2;
	double cs2 = sl*sl/3;

	for(int i=0;i<eq->size;i++){
		d2 = ux[i]*ux[i]+uy[i]*uy[i];
		for(int j=0;j<eq->nq;j++){
			idx = j+i*eq->nq;
			d1 = ux[i]*LC[0+2*j]+uy[i]*LC[1+2*j];
			eq->m[idx]=W[j]*D[i]*(1+d1/(cs2)+(d1*d1)/(cs2*cs2*2)-d2/(2*cs2));
		}
	}
	get_density(D2,eq);
	for(int i=0;i<eq->size;i++){
		for(int j=0;j<eq->nq;j++){
			idx = j+i*eq->nq;
			eq->m[idx] = eq->m[idx] * D[i] / D2[i];
		}
	}

	free(D);
	free(D2);
	free(ux);
	free(uy);
}

void ND_collision(struct ND *res, struct ND *old, struct ND *eq, struct BC *bc,  double cf){
	int idx_nd;
	int idx_nd_p;
	int idx_bc;
	for(int j=0;j<res->ny;j++){
		for(int i=0;i<res->nx;i++){
			idx_bc=i+1+(j+1)*bc->nx;
			for(int vc=0;vc<res->nq;vc++){
				idx_nd = vc+(i+j*res->nx)*res->nq;
				switch(bc->m[idx_bc]){
					case '+':
						res->m[idx_nd] = old->m[idx_nd]*(1-cf) + eq->m[idx_nd]*cf;
						break;
					case '>':
						res->m[idx_nd] = bc->s[vc];
						break;
					case '#':
						res->m[idx_nd] = old->m[idx_nd_p];
						break;
				}
			}
		}
	}
}

void ND_streaming(struct ND *res, struct ND *old, struct BC *bc, struct FC *fc){
	double d;
	int idx_nd, idx_nd_p, idx_nd_r;
	int idx_bc, idx_bc_p, idx_bc_px;

	regex_t not_wall;
	regcomp(&not_wall, "^[0-9]" ,0);

	for(int j=0;j<res->ny;j++){
		for(int i=0;i<res->nx;i++){
			idx_bc=i+1+(j+1)*bc->nx;

			for(int vc=0;vc<res->nq;vc++){
				idx_nd = vc+(i+j*res->nx)*res->nq;
				idx_nd_p = vc+((i - LC[0+vc*2])+(j - LC[1+vc*2])*res->nx)*res->nq;
				idx_nd_r = (8-vc) + (i+j*res->nx)*res->nq;
				idx_bc_p = (i - LC[0+vc*2] + 1) + (j - LC[1+vc*2] + 1) * bc->nx;

				if (bc->m[idx_bc] == '+' ){
					if(!regexec(&not_wall,&bc->m[idx_bc_p],0,NULL,0)){
						//transfer ascii number to int by minusing 48
						res->m[idx_nd] = old->m[idx_nd_r];
					}else{
						res->m[idx_nd] = old->m[idx_nd_p];
					}
				} else if (isnum(bc->m[idx_bc])){
					res->m[idx_nd] = 0;
				}else{
					res->m[idx_nd] = old->m[idx_nd];
				}
			}
		}
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

void jet_colormap(double num, double max, int *c){
	double gap = max/4;
	double slope = 125/gap;
	if( 0 <=num && num < gap){
		c[0] = (int)(125 - slope*num);
		c[1] = 250;
		c[2] = (int)(125 + slope*num);
	}else if ( gap <= num && num < 3*gap){
		c[0] = 0;
		c[1] = (int)(250 - slope*(num-gap));
		c[2] = 250;
	}else if ( 3*gap <= num){
		c[0] = 0;
		c[1] = 0;
		c[1] = (int)(250 - slope*(num-3*gap));
	}else if ( num < 0 && num >= -gap ){
		c[0] = (int)(125 + slope*-num);
		c[1] = 250;
		c[2] = (int)(125 - slope*-num);
	}else if ( num < -gap && num >= -3*gap ){
		c[0] = 250;
		c[1] = (int)(250 - slope*(-num-gap));
		c[2] = 0;
	}else if ( num < -3*gap){
		c[0] = (int)(250 - slope*(-num-3*gap));
		c[1] = 0;
		c[2] = 0;
	}
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
void nd_pgm_write(double *m, int nx, int ny, int pmax, int scale, FILE *f){
	fprintf(f,"P2\n");
	fprintf(f,"%d %d\n",nx*scale,ny*scale);
	fprintf(f,"255\n");
	for(int j=0;j<ny;j++){
		for(int sy=0;sy<scale;sy++){
			for(int i=0;i<nx;i++){
				for(int sx=0;sx<scale;sx++){
					fprintf(f,"%4d ",(int)(255*m[i+j*nx]/pmax));
				}
			}
			fprintf(f,"\n");
		}
	}
}

void nd_ppm_write(double *m, int nx, int ny, double pmax, int scale, FILE *f){
	fprintf(f,"P3\n");
	fprintf(f,"%d %d 255\n",nx*scale,ny*scale);
	int *v;
	for(int j=0;j<ny;j++){
		for(int sy=0;sy<scale;sy++){
			for(int i=0;i<nx;i++){
				for(int sx=0;sx<scale;sx++){
					v=(int *)malloc(3*sizeof(int));
					jet_colormap(m[i+j*nx],pmax,v);
					fprintf(f,"%d %d %d\n",v[0],v[1],v[2]);
					free(v);
				}
				int nx,ny;
			}
		}
	}
}

