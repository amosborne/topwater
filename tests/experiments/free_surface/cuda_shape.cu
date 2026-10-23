#include "cuda_dual.cuh"
#include "cuda_wide_dual.cuh"

#ifndef SHAPE_WIDTH
#define SHAPE_WIDTH 1
#endif
#if SHAPE_WIDTH == 1
using ShapeDual=Dual;
__device__ void set_lane(ShapeDual& a,int k,double value){a.d=value;}
__device__ double get_lane(ShapeDual a,int k){return a.d;}
#else
using ShapeDual=WideDual<SHAPE_WIDTH>;
__device__ void set_lane(ShapeDual& a,int k,double value){a.d[k]=value;}
__device__ double get_lane(ShapeDual a,int k){return a.d[k];}
#endif

struct DualView {
    const double *value,*jacobian;
    int stride,seed,side;
    __device__ ShapeDual operator[](int index) const {
        int face=index/stride, local=index%stride;
        ShapeDual result(value[index]);
        for(int k=0;k<SHAPE_WIDTH;++k){
            int lane=seed+k;
            set_lane(result,k,lane/9==side?jacobian[(face*9+lane%9)*stride+local]:0.);
        }
        return result;
    }
};

#include "cuda_pair.cuh"

extern "C" __global__ void vertex_gradient(
    const double *tri0,const double *triJ,const double *normal0,const double *normalJ,
    const double *outward0,const double *outwardJ,const double *tangent0,const double *tangentJ,
    const double *length0,const double *lengthJ,const double *area0,const double *areaJ,
    const double *center0,const double *centerJ,const double *sign0,const double *signJ,
    const double *gradient0,const double *gradientJ,const double *points0,const double *pointsJ,
    const double *weights,const long long *faces,const double *curlphi0,const double *curlphiJ,
    const double *phi,const double *aux,int nf,int nv,int nq,double density,double *out) {
    int work=blockDim.x*blockIdx.x+threadIdx.x;
    const int batches=18/SHAPE_WIDTH;
    if (work>=nf*nf*batches) return;
    int pair=work/batches,seed=(work%batches)*SHAPE_WIDTH,i=pair/nf,j=pair%nf;
    DualView tri={tri0,triJ,9,seed,1},normal={normal0,normalJ,3,seed,1};
    DualView outward={outward0,outwardJ,9,seed,1},tangent={tangent0,tangentJ,9,seed,1};
    DualView length={length0,lengthJ,3,seed,1},source_area={area0,areaJ,1,seed,1};
    DualView target_area={area0,areaJ,1,seed,0},center={center0,centerJ,3,seed,1};
    DualView sign={sign0,signJ,1,seed,1},gradient={gradient0,gradientJ,9,seed,1};
    DualView points={points0,pointsJ,nq*3,seed,0};
    DualView ci={curlphi0,curlphiJ,3,seed,0},cj={curlphi0,curlphiJ,3,seed,1};
    ShapeDual vsum,linear[3];
    energy_pair<ShapeDual,DualView>(tri,normal,outward,tangent,length,source_area,target_area,
        center,sign,gradient,points,weights,i,j,nq,vsum,linear);

    ShapeDual curl_product=0.;
    for(int k=0;k<3;++k) curl_product+=ci[i*3+k]*cj[j*3+k];
    ShapeDual scalar=(.5*density)*vsum*(curl_product-aux[i]*aux[j]);
    for(int c=0;c<3;++c)
        scalar+=density*aux[i]*linear[c]*phi[faces[j*3+c]];
    for(int k=0;k<SHAPE_WIDTH;++k){
        int lane=seed+k,local=lane%9,face=lane<9?i:j;
        atomicAdd(out+faces[face*3+local/3]*3+local%3,get_lane(scalar,k));
    }
}
