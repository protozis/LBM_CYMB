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
	1.414213,
	1,
	1.414213,
	1,
	0,
	1,
	1.414213,
	1,
	1.414213
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
__private double get_dist(uint *addr, __global double *b, uint nq){
	__private double tmp = 0;
	for(int i=0;i<nq;i++){
		tmp += pow((double)b[i] - (double)addr[i],2);
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
		dist = get_dist(addr,&bck[i*bc_nq*CY_KIE_NUM+4*bc_nq],bc_nq) - bcp[i*CY_PAR_NUM+3];
		if(dist <= 0){
			res = i;
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
	1.414213,
void get_intercept(uint *addr,uint vc, uint obj,uint bc_nq,__global double *bcp, __global double *bck, double *out){
}
__kernel void propagate(uint nq, double sl, double cf, __global double *nd, __global double *res, uint bc_no, uint bc_nq, __global double *bcv, __global double *bcp, __global double *bck){
	uint addr[] = {get_global_id(0),get_global_id(1)};
	uint addr_p[2];
	uint nx = get_global_size(0);
	uint ny = get_global_size(1);
	uint idx_nd,idx_nd_p;
	uint obj;
	double macro[3];
	double dist;
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
				idx_nd = (addr_p[0] + addr_p[1]*nx)*nq;
				obj = is_inside(addr_p,bc_no,bc_nq,bcp,bck);
				if(obj != NULL){
					dist = get_dist(addr,&bck[obj*bc_nq*CY_KIE_NUM+4*bc_nq],bc_nq);
					if(dist < LE[vc]/2){
					}else{

					}
				} else {
					res[vc+idx_nd] = nd[vc+idx_nd_p];
				}
			}
			get_macro(&res[idx_nd],nq,&macro[0]);
		}
	}
}
