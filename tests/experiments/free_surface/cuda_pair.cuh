// Single numerical SAMPLE formula for primal energy and shape derivatives.
// Adapters below only weight/contract these same integral moments.
struct NullSampleObserver {
    template<class S> __device__ void point(const S*,const S[3][3],const S*){}
    template<class S> __device__ void solid(S,S,S){}
    template<class S> __device__ void edge(int,S,S,S,S){}
};

template<class Scalar,class View,class Observer>
__device__ void energy_sample_observed(
    View tri,View normal,View outward,View tangent,View length,
    View source_area,View center,View sign,View points,int i,int j,int nq,int q,
    Scalar& single,Scalar& dbl,Scalar moment[3],Observer& observer) {
    const Scalar factor=4.*3.1415926535897932384626433832795;
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
        observer.point(x,r,distance);
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
            observer.solid(height,angle,denominator);
        }
        else observer.solid(height,angle,Scalar(0.));
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
            observer.edge(e,d,a,denominator,logarithm);
            terms+=d*logarithm;
            #pragma unroll
            for (int k=0;k<3;++k) edge_gradient[k]+=outward[j*9+e*3+k]*logarithm;
        }
        single=(terms-height*angle)/factor;
        dbl=angle*sign[j]/factor;
        for(int k=0;k<3;++k)moment[k]=0.;
        if (i!=j) {
            #pragma unroll
            for (int k=0;k<3;++k)
                moment[k]=(x[k]-center[j*3+k]-height*normal[j*3+k])*dbl
                    -height*edge_gradient[k]*sign[j]/factor;
        }
}

template<class Scalar,class View>
__device__ void energy_sample(
    View tri,View normal,View outward,View tangent,View length,
    View source_area,View center,View sign,View points,int i,int j,int nq,int q,
    Scalar& single,Scalar& dbl,Scalar moment[3]) {
    NullSampleObserver observer;
    energy_sample_observed<Scalar,View,NullSampleObserver>(tri,normal,outward,
        tangent,length,source_area,center,sign,points,i,j,nq,q,single,dbl,moment,observer);
}

// Existing raw V/H interface: unchanged quadrature, self branch and jump.
// Views use flat indices; source and target area may carry different seeds.
template<class Scalar,class View>
__device__ void energy_pair(
    View tri,View normal,View outward,View tangent,View length,
    View source_area,View target_area,View center,View sign,
    View gradient,View points,const double* weights,int i,int j,int nq,
    Scalar& vsum,Scalar linear[3]) {
    vsum=0.;for(int c=0;c<3;++c)linear[c]=0.;
    for(int q=0;q<nq;++q){
        Scalar single,dbl,moment[3];
        energy_sample<Scalar,View>(tri,normal,outward,tangent,length,source_area,
            center,sign,points,i,j,nq,q,single,dbl,moment);
        Scalar weight=target_area[i]*weights[q];
        vsum+=weight*single;
        if(i!=j){
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

// Same sample moments, contracted BEFORE constructing three local H columns.
// A and g remain shape-active; zi and nodal phi are stationary constants.
template<class Scalar,class View>
__device__ Scalar stationary_pair_scalar(
    View tri,View normal,View outward,View tangent,View length,
    View source_area,View target_area,View center,View sign,View points,
    const double* weights,int i,int j,int nq,Scalar A,const Scalar g[3],
    double zi,double sum_phi) {
    Scalar result=0.;
    for(int q=0;q<nq;++q){
        Scalar single,dbl,moment[3];
        energy_sample<Scalar,View>(tri,normal,outward,tangent,length,source_area,
            center,sign,points,i,j,nq,q,single,dbl,moment);
        Scalar value=A*single;
        if(i!=j){
            Scalar mixed=dbl*(sum_phi/3.);
            for(int k=0;k<3;++k)mixed+=moment[k]*g[k];
            value+=zi*mixed;
        }
        result+=target_area[i]*weights[q]*value;
    }
    if(i==j)result+=zi*target_area[i]*sum_phi/6.;
    return result;
}
