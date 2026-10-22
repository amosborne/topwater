extern "C" __global__ void energy_pairs(
    const double* tri, const double* normal, const double* outward,
    const double* tangent, const double* length, const double* area,
    const double* center, const double* sign, const long long* faces,
    const double* gradient, const double* points, const double* weights,
    int nf, int nv, int nq, double* v, double* h) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if (pair>=nf*nf) return;
    int i=pair/nf, j=pair%nf;
    double vsum=0., linear[3]={0.,0.,0.};
    const double factor=4.*3.1415926535897932384626433832795;
    for (int q=0;q<nq;++q) {
        double r[3][3], distance[3], edge_gradient[3]={0.,0.,0.};
        const double* x=points+(i*nq+q)*3;
        #pragma unroll
        for (int e=0;e<3;++e) {
            double squared=0.;
            #pragma unroll
            for (int k=0;k<3;++k) {
                r[e][k]=x[k]-tri[j*9+e*3+k];
                squared+=r[e][k]*r[e][k];
            }
            distance[e]=sqrt(squared);
        }
        double height=0.;
        #pragma unroll
        for (int k=0;k<3;++k) height+=r[0][k]*normal[j*3+k];
        if (i==j) height=0.;
        double angle=0.;
        if (i!=j) {
            double denominator=distance[0]*distance[1]*distance[2];
            #pragma unroll
            for (int e=0;e<3;++e) {
                double dot=0.;
                #pragma unroll
                for (int k=0;k<3;++k) dot+=r[e][k]*r[(e+1)%3][k];
                denominator+=dot*distance[(e+2)%3];
            }
            angle=2.*atan2(2.*area[j]*height,denominator);
        }
        double terms=0.;
        #pragma unroll
        for (int e=0;e<3;++e) {
            double d=0.,a=0.;
            #pragma unroll
            for (int k=0;k<3;++k) {
                d-=r[e][k]*outward[j*9+e*3+k];
                a+=r[e][k]*tangent[j*9+e*3+k];
            }
            double l=length[j*3+e],end=l-a,perp=d*d+height*height;
            double denominator=perp/(distance[e]+fabs(a))
                +perp/(distance[(e+1)%3]+fabs(end))
                +2.*(fmax(-a,0.)+fmax(a-l,0.));
            double logarithm=log1p(2.*l/denominator);
            terms+=d*logarithm;
            #pragma unroll
            for (int k=0;k<3;++k) edge_gradient[k]+=outward[j*9+e*3+k]*logarithm;
        }
        double single=(terms-height*angle)/factor;
        double dbl=angle*sign[j]/factor;
        double weight=area[i]*weights[q];
        vsum+=weight*single;
        if (i!=j) {
            double moment[3];
            #pragma unroll
            for (int k=0;k<3;++k)
                moment[k]=(x[k]-center[j*3+k]-height*normal[j*3+k])*dbl
                    -height*edge_gradient[k]*sign[j]/factor;
            #pragma unroll
            for (int c=0;c<3;++c) {
                double value=dbl/3.;
                #pragma unroll
                for (int k=0;k<3;++k) value+=moment[k]*gradient[j*9+c*3+k];
                linear[c]+=weight*value;
            }
        }
    }
    v[pair]=vsum;
    #pragma unroll
    for (int c=0;c<3;++c)
        atomicAdd(h+i*nv+faces[j*3+c],linear[c]+(i==j?area[i]/6.:0.));
}

extern "C" __global__ void curl_energy(
    const double* v, const double* curl, const int* indptr,
    const int* incident_face, const int* incident_corner,
    int nf, int nv, double* w) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if (pair>=nv*nv) return;
    int i=pair/nv,j=pair%nv;
    double sum=0.;
    for (int a=indptr[i];a<indptr[i+1];++a) {
        int fi=incident_face[a],ci=incident_corner[a];
        const double* left=curl+fi*9+ci*3;
        for (int b=indptr[j];b<indptr[j+1];++b) {
            int fj=incident_face[b],cj=incident_corner[b];
            const double* right=curl+fj*9+cj*3;
            sum+=v[fi*nf+fj]*(left[0]*right[0]+left[1]*right[1]+left[2]*right[2]);
        }
    }
    w[pair]=sum;
}
