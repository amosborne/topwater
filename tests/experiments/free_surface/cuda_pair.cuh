// Single numerical pair formula for primal energy and shape derivatives.
// Views use flat indices; source and target area views may carry different seeds.
template<class Scalar,class View>
__device__ void energy_pair(
    View tri,View normal,View outward,View tangent,View length,
    View source_area,View target_area,View center,View sign,
    View gradient,View points,const double* weights,int i,int j,int nq,
    Scalar& vsum,Scalar linear[3]) {
    vsum=0.; for (int c=0;c<3;++c) linear[c]=0.;
    const Scalar factor=4.*3.1415926535897932384626433832795;
    for (int q=0;q<nq;++q) {
        Scalar r[3][3], distance[3], edge_gradient[3]={0.,0.,0.};
        Scalar x[3]; for (int k=0;k<3;++k) x[k]=points[(i*nq+q)*3+k];
        #pragma unroll
        for (int e=0;e<3;++e) {
            Scalar squared=0.;
            #pragma unroll
            for (int k=0;k<3;++k) {
                r[e][k]=x[k]-tri[j*9+e*3+k];
                squared+=r[e][k]*r[e][k];
            }
            distance[e]=sqrt(squared);
        }
        Scalar height=0.;
        #pragma unroll
        for (int k=0;k<3;++k) height+=r[0][k]*normal[j*3+k];
        if (i==j) height=0.;
        Scalar angle=0.;
        if (i!=j) {
            Scalar denominator=distance[0]*distance[1]*distance[2];
            #pragma unroll
            for (int e=0;e<3;++e) {
                Scalar dot=0.;
                #pragma unroll
                for (int k=0;k<3;++k) dot+=r[e][k]*r[(e+1)%3][k];
                denominator+=dot*distance[(e+2)%3];
            }
            angle=2.*atan2(2.*source_area[j]*height,denominator);
        }
        Scalar terms=0.;
        #pragma unroll
        for (int e=0;e<3;++e) {
            Scalar d=0.,a=0.;
            #pragma unroll
            for (int k=0;k<3;++k) {
                d-=r[e][k]*outward[j*9+e*3+k];
                a+=r[e][k]*tangent[j*9+e*3+k];
            }
            Scalar l=length[j*3+e],end=l-a,perp=d*d+height*height;
            Scalar denominator=perp/(distance[e]+fabs(a))
                +perp/(distance[(e+1)%3]+fabs(end))
                +2.*(fmax(-a,0.)+fmax(a-l,0.));
            Scalar logarithm=log1p(2.*l/denominator);
            terms+=d*logarithm;
            #pragma unroll
            for (int k=0;k<3;++k) edge_gradient[k]+=outward[j*9+e*3+k]*logarithm;
        }
        Scalar single=(terms-height*angle)/factor;
        Scalar dbl=angle*sign[j]/factor;
        Scalar weight=target_area[i]*weights[q];
        vsum+=weight*single;
        if (i!=j) {
            Scalar moment[3];
            #pragma unroll
            for (int k=0;k<3;++k)
                moment[k]=(x[k]-center[j*3+k]-height*normal[j*3+k])*dbl
                    -height*edge_gradient[k]*sign[j]/factor;
            #pragma unroll
            for (int c=0;c<3;++c) {
                Scalar value=dbl/3.;
                #pragma unroll
                for (int k=0;k<3;++k) value+=moment[k]*gradient[j*9+c*3+k];
                linear[c]+=weight*value;
            }
        }
    }
    // Smooth self double-layer moments vanish; the mixed trace jump occurs once.
    for (int c=0;c<3;++c) linear[c]+=(i==j?target_area[i]/6.:Scalar(0.));
}
