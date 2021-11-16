#define FC_OFFSET 1000
#define CS_LTTC_2 0.333333
#define CS_LTTC 0.57735
#define FCCS 5773.5
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable
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
__constant double LV[18] = {
	-0.707106,-0.707106,
	0,-1,
	0.707106,-0.707106,
	-1,0,
	0,0,
	1,0,
	-0.707106,0.707106,
	0,1,
	0.707106,0.707106,
};
__constant double LE[9] = {
	1.414214,
	1,
	1.414214,
	1,
	0,
	1,
	1.414214,
	1,
	1.414214
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
__constant double WA[9] = {
	0.0555556,
	0.2222222,
	0.0555556,
	0.2222222,
	0,
	0.2222222,
	0.0555556,
	0.2222222,
	0.0555556
};
double get_dist_int(int *addr, double *b, int nq){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += pow((b[i] - (double)addr[i]),2);
	}
	return sqrt(tmp);
}
double get_dist_int_global(int *addr, __global double *b, int nq){
	double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += pow((b[i] - (double)addr[i]),2);
	}
	return sqrt(tmp);
}
int is_border(int *addr, int nx, int ny){
	return addr[0] == 0 || addr[0] == nx-1 || addr[1] == 0 || addr[1] == ny-1;
}
int is_inside(int *addr, int bc_no, int bc_nq, __global double *bcpos,__global double *bcrad){
	double dist;
	int res = 0;
	for(int i=0;i<bc_no;i++){
		dist = get_dist_int_global(addr,&bcpos[i*bc_nq],bc_nq) - bcrad[i];
		if(dist <= 0){
			res = i+1;
			break;
		}
	}
	return res;
}
void get_macro(__global double *nd, int nq, __private double *out){
	out[0] = 0;
	out[1] = 0;
	out[2] = 0;
	for(int i=0;i<nq;i++){
		out[0] += nd[i];
		out[1] += nd[i] * LC[0+2*i];
		out[2] += nd[i] * LC[1+2*i];
	}
	if(out[0] != 0){
		out[1] = out[1]/out[0];
		out[2] = out[2]/out[0];
	} else {
		out[1] = 0;
		out[2] = 0;
	}
}
void solver(double a, double b, double c, double *res){
	res[0] = (-1*b + sqrt(b*b - 4*a*c))/(2*a);
	res[1] = (-1*b - sqrt(b*b - 4*a*c))/(2*a);
}
double border_dist(int *addr,int vc,int obj,int bc_nq, __global double *bcpos, __global double *bcrad){
	double r = bcrad[obj];
	double addrd[2] = {
		(double)addr[0],
		(double)addr[1]
	};
	double cir[2];
	cir[0] = bcpos[obj*bc_nq];
	cir[1] = bcpos[1+obj*bc_nq];
	double res[2];
	double icp[2];
	double a,b,c,d;
	switch(vc){
		case 0:
			d = addrd[1] - addrd[0] - cir[1];
			a = 2;
			b = 2*(d-cir[0]);
			c = pow(cir[0],2) + pow(d,2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[0];
			icp[1] = res[0] + addrd[1] - addrd[0];
			break;
		case 8:
			d = addrd[1] - addrd[0] - cir[1];
			a = 2;
			b = 2*(d-cir[0]);
			c = pow(cir[0],2) + pow(d,2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[1];
			icp[1] = res[1] + addrd[1] - addrd[0];
			break;
		case 1:
			a = 1;
			b = -2*cir[1];
			c = pow(cir[1],2) + pow((addrd[0]-cir[0]),2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = addrd[0];
			icp[1] = res[0];
			break;
		case 7:
			a = 1;
			b = -2*cir[1];
			c = pow(cir[1],2) + pow((addrd[0]-cir[0]),2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = addrd[0];
			icp[1] = res[1];
			break;
		case 2:
			d = addrd[1] + addrd[0] - cir[1];
			a = 2;
			b = -2*(d+cir[0]);
			c = pow(cir[0],2) + pow(d,2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[1];
			icp[1] = addrd[1] + addrd[0] - icp[0];
			break;
		case 6:
			d = addrd[1] + addrd[0] - cir[1];
			a = 2;
			b = -2*(d+cir[0]);
			c = pow(cir[0],2) + pow(d,2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[0];
			icp[1] = addrd[1] + addrd[0] - icp[0];
			break;
		case 3:
			a = 1;
			b = -2*cir[0];
			c = pow(cir[0],2) + pow((addrd[1]-cir[1]),2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[0];
			icp[1] = addrd[1];
			break;
		case 5:
			a = 1;
			b = -2*cir[0];
			c = pow(cir[0],2) + pow((addrd[1]-cir[1]),2) - pow(r,2);
			solver(a,b,c,res);
			icp[0] = res[1];
			icp[1] = addrd[1];
			break;
	}
	//printf("%lf %lf\n",icp[0],icp[1]);
	return get_dist_int(addr,icp,bc_nq)/LE[vc];
}
double dot_product(__constant double *a,__global double *b,int nq){
	double tmp;
	for(int i=0;i<nq;i++){
		tmp += a[i]*b[i];
	}
	return tmp;
}
void get_eq(__private double *eq, int nq, double *mo){
	double v1, v2;
	v2 = mo[1]*mo[1] + mo[2]*mo[2];
	for (int vc=0;vc<nq;vc++){
		v1 = mo[1]*LC[0+2*vc] + mo[2]*LC[1+2*vc];
		eq[vc] = W[vc]*mo[0]*(1 + v1/(CS_LTTC_2) + (v1*v1)/(CS_LTTC_2*CS_LTTC_2*2) - v2/(2*CS_LTTC_2));
	}
}
void addr_vec(int *ori, int *dst, int vc, int q){
	for(int i=0;i<2;i++){
		dst[i] = ori[i] + q*(int)LC[i+vc*2];
	}
}
int addr_idx(int *addr,int nx, int nq){
	return (addr[0] + addr[1]*nx)*nq;
}
__kernel void propagate(int nq, double cf, __global double *nd, __global double *res, int bc_no, int bc_nq, __global double *bcv, __global double *bcpos,__global double *bcpos_p, __global double *bcvel, __global double *bcrad, __global long *bcfc, double refuel_rto, double eat_rto){
	int addr[] = {get_global_id(0),get_global_id(1)};
	int addr_p[2];
	int addr_f[2];
	int addr_eat[2];
	int nx = get_global_size(0);
	int ny = get_global_size(1);
	int idx_nd,idx_nd_p,idx_nd_f;
	int obj_vc_p;
	double dq;
	int spc[9];
	__private double macro[3];
	double eq[9];
	idx_nd = addr_idx(addr,nx,nq);
	if(is_border(addr,nx,ny)){
		for(int vc=0;vc<nq;vc++){
			res[vc+idx_nd] = bcv[vc];
		}
	} else {
		if(is_inside(addr,bc_no,bc_nq,bcpos,bcrad)){
			for(int vc=0;vc<nq;vc++){
				res[vc+idx_nd] = 0;
			}
		}else{
			for(int vc=0;vc<nq;vc++){
				addr_vec(addr,addr_p,8-vc,1);
				obj_vc_p = is_inside(addr_p,bc_no,bc_nq,bcpos,bcrad);
				idx_nd_p = addr_idx(addr_p,nx,nq);
				if(obj_vc_p != 0){
					if(!is_inside(addr_p,bc_no,bc_nq,bcpos_p,bcrad)){
						//eat
						res[vc+idx_nd] = nd[8-vc+idx_nd]+nd[vc+idx_nd_p]*eat_rto;
					}else if(is_inside(addr,bc_no,bc_nq,bcpos_p,bcrad)){
						//refuel
						addr_vec(addr,addr_f,vc,1);
						idx_nd_f = addr_idx(addr_f,nx,nq);
						if(is_inside(addr_p,bc_no,bc_nq,bcpos_p,bcrad)){
							res[vc+idx_nd] = (nd[8-vc+idx_nd_f]+nd[vc+idx_nd_f])*refuel_rto;
						}else{
							res[vc+idx_nd] = nd[vc+idx_nd_p]*refuel_rto;
						}
					}else{
						res[vc+idx_nd] = nd[8-vc+idx_nd];
					}
					for(int i=0;i<bc_nq;i++){
						atomic_add(&bcfc[i+(obj_vc_p-1)*bc_no],(int)(-2*(res[vc+idx_nd])*LV[i+vc*2]*LE[vc]*FCCS));
					}
				} else {
					res[vc+idx_nd] = nd[vc+idx_nd_p];
				}
			}
			get_macro(&res[idx_nd],nq,macro);
			get_eq(eq,nq,macro);
			for(int vc=0;vc<nq;vc++){
				res[vc+idx_nd] = res[vc+idx_nd]*(1-cf) + eq[vc]*cf;
			}
		}
	}
}
