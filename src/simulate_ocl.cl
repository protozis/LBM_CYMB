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
__private double get_dist_uint(uint *addr, double *b, uint nq){
	__private double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += pow(((double)b[i] - (double)addr[i]),2);
	}
	return sqrt(tmp);
}
uint is_border(uint *addr, uint nx, uint ny){
	return addr[0] == 0 || addr[0] == nx-1 || addr[1] == 0 || addr[1] == ny-1;
}
uint is_inside(uint *addr, uint bc_no, uint bc_nq, __global double *bcp,__global double *bck){
	double dist;
	uint res = NULL;
	for(int i=0;i<bc_no;i++){
		dist = get_dist_uint(addr,&bck[i*bc_nq*CY_KIE_NUM+4*bc_nq],bc_nq) - bcp[i*CY_PAR_NUM+3];
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
double border_dist(uint *addr,uint vc,uint obj,uint bc_nq, __global double *bcp, __global double *bck){
	double r = bcp[obj*CY_PAR_NUM+3];
	double addrd[2] = {
		(double)addr[0],
		(double)addr[1]
	};
	double *cir = &bck[obj*bc_nq*CY_KIE_NUM+4*bc_nq];
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
//	printf("%lf %lf %lf %lf %lf %lf\n",addrd[0],addrd[1],LC[vc*2],LC[1+vc*2],icp[0],icp[1]);
	return get_dist_uint(addr,icp,bc_nq)/LE[vc];
}
double dot_product(__constant double *a,double *b,uint nq){
	double tmp;
	for(int i=0;i<nq;i++){
		tmp += a[i]*b[i];
	}
	return tmp;
}
void get_eq(__private double *eq, int nq, double sl, double *mo){
	double v1, v2;
	double cs2 = sl*sl/3;
	for (int vc=0;vc<nq;vc++){
		v1 = mo[1]*LC[0+2*vc] + mo[2]*LC[1+2*vc];
		v2 = mo[1]*mo[1] + mo[2]*mo[2];
		eq[vc] = W[vc]*mo[0]*(1 + v1/(cs2) + (v1*v1)/(cs2*cs2*2) - v2/(2*cs2));
	}
}
__kernel void propagate(uint nq, double sl, double cf, __global double *nd, __global double *res, uint bc_no, uint bc_nq, __global double *bcv, __global double *bcp, __global double *bck){
	uint addr[] = {get_global_id(0),get_global_id(1)};
	uint addr_p[2];
	uint nx = get_global_size(0);
	uint ny = get_global_size(1);
	uint idx_nd,idx_nd_p,idx_nd_f,idx_nd_ff;
	uint objp = NULL;
	double dq;
	double macro[3];
	double eq[9];
	idx_nd = (addr[0] + addr[1]*nx)*nq;
	if(is_border(addr,nx,ny)){
		for(int vc=0;vc<nq;vc++){
			res[vc+idx_nd] = bcv[vc];
		}
	} else {
		if(is_inside(addr,bc_no,bc_nq,bcp,bck)){
			for(int vc=0;vc<nq;vc++){
				res[vc+idx_nd] = 0;
			}
		}else{
			for(int vc=0;vc<nq;vc++){
				addr_p[0] = addr[0] - LC[0+vc*2]; 
				addr_p[1] = addr[1] - LC[1+vc*2];
				objp = is_inside(addr_p,bc_no,bc_nq,bcp,bck);
				if(objp != 0){
					if(vc != 4){
						dq = border_dist(addr,8-vc,objp-1,bc_nq,bcp,bck);
						idx_nd_f = ((addr[0] + LC[0+vc*2]) + (addr[1]+LC[1+vc*2])*nx)*nq;
						idx_nd_ff = ((addr[0] + 2*LC[0+vc*2]) + (addr[1]+2*LC[1+vc*2])*nx)*nq;
						if(dq < 0.5){
							res[vc+idx_nd] = dq*(1+2*dq)*nd[8-vc+idx_nd] + (1 - 4*dq*dq)*nd[8-vc+idx_nd_f] - dq*(1-2*dq)*nd[8-vc+idx_nd_ff] + 3*WA[8-vc]*dot_product(&LC[8-vc*2],&bck[(objp-1)*bc_nq*CY_KIE_NUM+2*bc_nq],bc_nq);
						}else{
							res[vc+idx_nd] = nd[8-vc+idx_nd]/(dq*(1+2*dq)) + (2*dq-1)*nd[8-vc+idx_nd_f]/dq - (2*dq-1)*nd[8-vc+idx_nd_ff]/(2*dq+1) + 3*WA[8-vc]*dot_product(&LC[8-vc*2],&bck[(objp-1)*bc_nq*CY_KIE_NUM+2*bc_nq],bc_nq)/(dq*(2*dq+1));
						}
					}
				} else {
					idx_nd_p = (addr_p[0] + addr_p[1]*nx)*nq;
					res[vc+idx_nd] = nd[vc+idx_nd_p];
				}
			}
			get_macro(&res[idx_nd],nq,macro);
			get_eq(eq,nq,sl,macro);
			for(int vc=0;vc<nq;vc++){
				res[vc+idx_nd] = res[vc+idx_nd]*(1-cf) + eq[vc]*cf;
			}
		}
	}
}
