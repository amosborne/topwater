struct DirectionalView {
    const double *value,*derivative;
    __device__ Dual operator[](int index) const {
        return Dual(value[index],derivative[index]);
    }
};

extern "C" __global__ void directional_pairs(
    const double *tri0,const double *triD,const double *normal0,const double *normalD,
    const double *outward0,const double *outwardD,const double *tangent0,const double *tangentD,
    const double *length0,const double *lengthD,const double *area0,const double *areaD,
    const double *center0,const double *centerD,const double *sign0,const double *signD,
    const double *gradient0,const double *gradientD,const double *points0,const double *pointsD,
    const double *weights,const long long *faces,int nf,int nv,int nq,double *dv,double *dh) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if (pair>=nf*nf) return;
    int i=pair/nf,j=pair%nf;
    DirectionalView tri={tri0,triD},normal={normal0,normalD};
    DirectionalView outward={outward0,outwardD},tangent={tangent0,tangentD};
    DirectionalView length={length0,lengthD},area={area0,areaD};
    DirectionalView center={center0,centerD},sign={sign0,signD};
    DirectionalView gradient={gradient0,gradientD},points={points0,pointsD};
    Dual single,linear[3];
    energy_pair<Dual,DirectionalView>(tri,normal,outward,tangent,length,area,area,
        center,sign,gradient,points,weights,i,j,nq,single,linear);
    dv[pair]=single.d;
    for (int c=0;c<3;++c) atomicAdd(dh+i*nv+faces[j*3+c],linear[c].d);
}
