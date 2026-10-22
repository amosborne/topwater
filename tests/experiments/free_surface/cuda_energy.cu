extern "C" __global__ void energy_pairs(
    const double* tri, const double* normal, const double* outward,
    const double* tangent, const double* length, const double* area,
    const double* center, const double* sign, const long long* faces,
    const double* gradient, const double* points, const double* weights,
    int nf, int nv, int nq, double* v, double* h) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if (pair>=nf*nf) return;
    int i=pair/nf, j=pair%nf;
    double vsum,linear[3];
    energy_pair<double,const double*>(tri,normal,outward,tangent,length,area,area,
        center,sign,gradient,points,weights,i,j,nq,vsum,linear);
    v[pair]=vsum;
    for (int c=0;c<3;++c) atomicAdd(h+i*nv+faces[j*3+c],linear[c]);
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
