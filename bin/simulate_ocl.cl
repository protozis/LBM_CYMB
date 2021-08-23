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

void atomic_add_2dFC( __global int *sig, uint vc, uint nq, uint ofs, double val){
	for(int i=0;i<nq;i++){
		atom_add(&sig[i], (int)(LC[i+vc*2] * val * ofs));
	}
}


__kernel void propagate(uint nq, double sl, double cf, __global double *nd, __global double *res, __global char *bc, __global double *bcv, uint fc_no, uint fc_nq, uint fc_ofs, __global int *force){

	uint idx_nd, idx_nd_p, idx_nd_r;
	uint idx_bc, idx_bc_p;
	double d,deq,ux,uy;

	uint addr_x = get_global_id(0);
	uint addr_y = get_global_id(1);
	uint nx = get_global_size(0);
	uint ny = get_global_size(1);
	uint nbx = nx+2;
	uint nby = ny+2;
	__private double eq[9];

	idx_bc=addr_x+1+(addr_y+1)*nbx;
	switch(bc[idx_bc]){
		case '+':
			for(int vc=0;vc<nq;vc++){
				idx_nd = vc + (addr_x + addr_y*nx)*nq;
				idx_nd_p = vc+((addr_x - LC[0+vc*2])+(addr_y - LC[1+vc*2])*nx)*nq;
				idx_nd_r = (8-vc) + (addr_x + addr_y*nx)*nq;
				idx_bc_p = (addr_x - LC[0+vc*2] + 1) + (addr_y - LC[1+vc*2] + 1) * nbx;

				if(isnum(bc[idx_bc_p])){
					atomic_add_2dFC(&force[0+fc_nq*(bc[idx_bc_p]-48)],vc,fc_nq,fc_ofs, -2 * nd[idx_nd_r]);
					res[idx_nd] = nd[idx_nd_r];
				}else{
					res[idx_nd] = nd[idx_nd_p];
				}
			}
			d = get_d_global(&res[(addr_x + addr_y*nx)*nq], nq);
			ux = get_ux(&res[(addr_x + addr_y*nx)*nq], nq, d);
			uy = get_uy(&res[(addr_x + addr_y*nx)*nq], nq, d);
			get_eq(eq, nq, sl, d, ux, uy);

			deq = get_d_local(&eq[0],nq);
			for(int vc=0;vc<nq;vc++){
				idx_nd = vc + (addr_x + addr_y*nx)*nq;
				res[idx_nd] = res[idx_nd]*(1-cf) + eq[vc]*d*cf/deq;
				nd[idx_nd] = res[idx_nd];
			}
			break;
		case '>':
			for(int vc=0;vc<nq;vc++){
				idx_nd = vc + (addr_x + addr_y*nx)*nq;
				res[idx_nd] = bcv[vc];
				nd[idx_nd] = res[idx_nd];
			}
			break;
		case '#':
			for(int vc=0;vc<nq;vc++){
				idx_nd = vc + (addr_x + addr_y*nx)*nq;
				idx_nd_p = vc+((addr_x - LC[0+vc*2])+(addr_y - LC[1+vc*2])*nx)*nq;
				res[idx_nd] = nd[idx_nd_p];
				nd[idx_nd] = res[idx_nd];
			}
			break;
	}
}	
