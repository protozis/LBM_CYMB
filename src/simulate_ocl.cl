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


__kernel void propagate(uint nq, double sl, double cf, __global double *nd, __global double *res, __global double *bcp, __global double *bck){
}	
