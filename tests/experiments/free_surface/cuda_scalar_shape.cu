// Compiled after the existing shape geometry view and shared pair source.
extern "C" __global__ void contracted_vertex_gradient(
    const double *tri0,const double *triJ,const double *normal0,const double *normalJ,
    const double *outward0,const double *outwardJ,const double *tangent0,const double *tangentJ,
    const double *length0,const double *lengthJ,const double *area0,const double *areaJ,
    const double *center0,const double *centerJ,const double *sign0,const double *signJ,
    const double *gradient0,const double *gradientJ,const double *points0,const double *pointsJ,
    const double *weights,const long long *faces,const double *curlphi0,const double *curlphiJ,
    const double *phi,const double *aux,int nf,int nv,int nq,double density,double *out
#ifdef PAIR_DIAGNOSTIC
    ,const long long *selected,int selected_count
#endif
    ) {
    int work=blockDim.x*blockIdx.x+threadIdx.x;
    const int batches=18/SHAPE_WIDTH;
#ifdef PAIR_DIAGNOSTIC
    if(work>=selected_count*batches)return;
    int slot=work/batches,pair=selected[slot];
#else
    if(work>=nf*nf*batches)return;
    int pair=work/batches;
#endif
    int seed=(work%batches)*SHAPE_WIDTH,i=pair/nf,j=pair%nf;
    DualView tri={tri0,triJ,9,seed,1},normal={normal0,normalJ,3,seed,1};
    DualView outward={outward0,outwardJ,9,seed,1},tangent={tangent0,tangentJ,9,seed,1};
    DualView length={length0,lengthJ,3,seed,1},source_area={area0,areaJ,1,seed,1};
    DualView target_area={area0,areaJ,1,seed,0},center={center0,centerJ,3,seed,1};
    DualView sign={sign0,signJ,1,seed,1},gradient={gradient0,gradientJ,9,seed,1};
    DualView points={points0,pointsJ,nq*3,seed,0};
    DualView ci={curlphi0,curlphiJ,3,seed,0},cj={curlphi0,curlphiJ,3,seed,1};
    ShapeDual A=-aux[i]*aux[j],g[3]={0.,0.,0.};double sum_phi=0.;
    for(int k=0;k<3;++k)A+=ci[i*3+k]*cj[j*3+k];
    A=A*.5;
    for(int c=0;c<3;++c){
        double value=phi[faces[j*3+c]];sum_phi+=value;
        for(int k=0;k<3;++k)g[k]+=value*gradient[j*9+c*3+k];
    }
    ShapeDual scalar=density*stationary_pair_scalar<ShapeDual,DualView>(tri,normal,
        outward,tangent,length,source_area,target_area,center,sign,points,weights,
        i,j,nq,A,g,aux[i],sum_phi);
    for(int k=0;k<SHAPE_WIDTH;++k){
        int lane=seed+k,local=lane%9,face=lane<9?i:j;
#ifdef PAIR_DIAGNOSTIC
        out[slot*19+lane]=get_lane(scalar,k);
#else
        atomicAdd(out+faces[face*3+local/3]*3+local%3,get_lane(scalar,k));
#endif
    }
#ifdef PAIR_DIAGNOSTIC
    if(seed==0)out[slot*19+18]=scalar.v;
#endif
}
