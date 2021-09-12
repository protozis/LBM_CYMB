#define CY_PAR_NUM 4
#define CY_KIE_NUM 5

__constant double LC[18] = {
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

__constant double W[9] ={ 
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
int isnum(char c){
	int tmp = 0;
	if( c >= '0' && c <= '9' ){
		tmp = 1;
	}
	return tmp;
}
double get_d_global(__global double *nd, int nq){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += nd[i];
	}
	return tmp;
}

double get_d_local(__private double *nd, int nq){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += nd[i];
	}
	return tmp;
}

double get_ux(__global double *nd, int nq , double d){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += nd[i] * LC[0+2*i];
	}
	if(d != 0){
		tmp = tmp/d;
	}else{
		tmp = 0;
	}
	return tmp;
}

double get_uy(__global double *nd, int nq , double d){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += nd[i] * LC[1+2*i];
	}
	if(d != 0){
		tmp = tmp/d;
	}else{
		tmp = 0;
	}
	return tmp;
}

void get_eq(__private double *eq, int nq, double sl, double d, double ux, double uy){
	double v1, v2;
	double cs2 = sl*sl/3;
	for (int vc=0;vc<nq;vc++){
		v1 = ux*LC[0+2*vc] + uy*LC[1+2*vc];
		v2 = ux*ux + uy*uy;
		eq[vc] = W[vc]*d*(1 + v1/(cs2) + (v1*v1)/(cs2*cs2*2) - v2/(2*cs2));
	}
}
__private double get_dist(__private double *a, __private double *b, uint nq){
	__private double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += pow(b[i] - a[i],2);
	}
	return sqrt(tmp);
}
__kernel void propagate(uint nq, double sl, double cf, __global double *nd, __global double *res, uint bc_no, uint bc_nq, __global double *bcv, __global double *bcp, __global double *bck){
	uint addr[] = {get_global_id(0),get_global_id(1)};
	uint nx = get_global_size(0);
	uint ny = get_global_size(1);
	uint idx_bck;

	double *bck_force, *bck_acc, *bck_vel, *bck_dsp, *bck_pos;
	for(int i=0;i<bc_no;i++){
		idx_bck = i*bc_nq*CY_KIE_NUM;
		bck_force = &bck[idx_bck+0*bc_nq];
		bck_acc = &bck[idx_bck+1*bc_nq];
		bck_vel = &bck[idx_bck+2*bc_nq];
		bck_dsp = &bck[idx_bck+3*bc_nq];
		bck_pos = &bck[idx_bck+4*bc_nq];
		printf("%lf %lf\n",bck_pos[0],bck_pos[1]);
	}
	
}	
