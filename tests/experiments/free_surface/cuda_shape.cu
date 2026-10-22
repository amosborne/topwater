#include "cuda_dual.cuh"

struct DualView {
    const double *value,*jacobian;
    int stride,seed;
    __device__ Dual operator[](int index) const {
        int face=index/stride, local=index%stride;
        return Dual(value[index],seed<0?0.:jacobian[(face*9+seed)*stride+local]);
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
    if (work>=nf*nf*18) return;
    int pair=work/18,lane=work%18,i=pair/nf,j=pair%nf;
    int target_seed=lane<9?lane:-1,source_seed=lane>=9?lane-9:-1;
    DualView tri={tri0,triJ,9,source_seed},normal={normal0,normalJ,3,source_seed};
    DualView outward={outward0,outwardJ,9,source_seed},tangent={tangent0,tangentJ,9,source_seed};
    DualView length={length0,lengthJ,3,source_seed},source_area={area0,areaJ,1,source_seed};
    DualView target_area={area0,areaJ,1,target_seed},center={center0,centerJ,3,source_seed};
    DualView sign={sign0,signJ,1,source_seed},gradient={gradient0,gradientJ,9,source_seed};
    DualView points={points0,pointsJ,nq*3,target_seed};
    DualView ci={curlphi0,curlphiJ,3,target_seed},cj={curlphi0,curlphiJ,3,source_seed};
    Dual vsum,linear[3];
    energy_pair<Dual,DualView>(tri,normal,outward,tangent,length,source_area,target_area,
        center,sign,gradient,points,weights,i,j,nq,vsum,linear);

    Dual curl_product=0.;
    for(int k=0;k<3;++k) curl_product+=ci[i*3+k]*cj[j*3+k];
    Dual scalar=(.5*density)*vsum*(curl_product-aux[i]*aux[j]);
    for(int c=0;c<3;++c)
        scalar+=density*aux[i]*linear[c]*phi[faces[j*3+c]];
    int seed=lane%9,face=lane<9?i:j;
    atomicAdd(out+faces[face*3+seed/3]*3+seed%3,scalar.d);
}
