// Prepended with the existing cuda_shape.cu definitions through energy_pair.
extern "C" __global__ void anchor_pairs(
    const double *tri0,const double *triJ,const double *normal0,const double *normalJ,
    const double *outward0,const double *outwardJ,const double *tangent0,const double *tangentJ,
    const double *length0,const double *lengthJ,const double *area0,const double *areaJ,
    const double *center0,const double *centerJ,const double *sign0,const double *signJ,
    const double *gradient0,const double *gradientJ,const double *points0,const double *pointsJ,
    const double *weights,int nf,int nq,double *value,double *jacobian) {
    int work=blockDim.x*blockIdx.x+threadIdx.x;
    const int batches=18/SHAPE_WIDTH;
    if(work>=nf*nf*batches)return;
    int pair=work/batches,seed=(work%batches)*SHAPE_WIDTH,i=pair/nf,j=pair%nf;
    DualView tri={tri0,triJ,9,seed,1},normal={normal0,normalJ,3,seed,1};
    DualView outward={outward0,outwardJ,9,seed,1},tangent={tangent0,tangentJ,9,seed,1};
    DualView length={length0,lengthJ,3,seed,1},source_area={area0,areaJ,1,seed,1};
    DualView target_area={area0,areaJ,1,seed,0},center={center0,centerJ,3,seed,1};
    DualView sign={sign0,signJ,1,seed,1},gradient={gradient0,gradientJ,9,seed,1};
    DualView points={points0,pointsJ,nq*3,seed,0};
    ShapeDual v,h[3];
    energy_pair<ShapeDual,DualView>(tri,normal,outward,tangent,length,source_area,
        target_area,center,sign,gradient,points,weights,i,j,nq,v,h);
    ShapeDual raw[4]={v,h[0],h[1],h[2]};
    for(int a=0;a<4;++a){
        if(seed==0)value[pair*4+a]=raw[a].v;
        for(int k=0;k<SHAPE_WIDTH;++k)
            jacobian[(pair*18+seed+k)*4+a]=get_lane(raw[a],k);
    }
}

extern "C" __global__ void linear_pairs(const double *value,const double *jacobian,
    const double *delta,const long long *faces,int nf,int nv,double *v,double *h) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if(pair>=nf*nf)return;
    int i=pair/nf,j=pair%nf;
    double raw[4];for(int a=0;a<4;++a)raw[a]=value[pair*4+a];
    for(int lane=0;lane<18;++lane){
        int face=lane<9?i:j,local=lane%9;
        double dx=delta[faces[face*3+local/3]*3+local%3];
        for(int a=0;a<4;++a)raw[a]+=jacobian[(pair*18+lane)*4+a]*dx;
    }
    v[pair]=raw[0];
    for(int c=0;c<3;++c)atomicAdd(h+i*nv+faces[j*3+c],raw[c+1]);
}

extern "C" __global__ void cached_pair_gradient(const double *jacobian,
    const long long *faces,const double *curlphi,const double *phi,const double *q,
    int nf,double density,double *out) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if(pair>=nf*nf)return;
    int i=pair/nf,j=pair%nf;
    double product=-q[i]*q[j];
    for(int a=0;a<3;++a)product+=curlphi[i*3+a]*curlphi[j*3+a];
    // product is symmetric: contraction with raw dV equals that with sym(dV).
    for(int lane=0;lane<18;++lane){
        const double *d=jacobian+(pair*18+lane)*4;
        double scalar=.5*product*d[0];
        for(int c=0;c<3;++c)scalar+=q[i]*d[c+1]*phi[faces[j*3+c]];
        int face=lane<9?i:j,local=lane%9;
        atomicAdd(out+faces[face*3+local/3]*3+local%3,density*scalar);
    }
}

extern "C" __global__ void current_curl_gradient(const long long *faces,
    const double *dc,const double *vc,int nf,double density,double *out) {
    int work=blockDim.x*blockIdx.x+threadIdx.x;
    if(work>=nf*9)return;
    int face=work/9,local=work%9;
    double scalar=0.;
    for(int a=0;a<3;++a)scalar+=dc[work*3+a]*vc[face*3+a];
    atomicAdd(out+faces[face*3+local/3]*3+local%3,density*scalar);
}
