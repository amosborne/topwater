// Uses the SAME analytic sample primitive as energy and shape derivatives.
// This entry point has explicit nonself ownership for arbitrary graph targets.
struct SourceGuard {
    double diameter, area, height, angle, angle_error, absolute_terms;
    bool inside;
    unsigned char flags;
    SourceGuard() = delete;
    __device__ SourceGuard(double d,double a):diameter(d),area(a),height(0.),
        angle(0.),angle_error(0.),absolute_terms(0.),inside(true),flags(0){}
    __device__ void point(const double*,const double r[3][3],const double* distance) {
        for(int e=0;e<3;++e) if(distance[e]==0.) flags|=1;
        double p[3]={0.,0.,0.};
        for(int e=0;e<3;++e) for(int k=0;k<3;++k) p[e]+=r[e][k]*r[(e+1)%3][k];
        angle_error=32.*2.2204460492503131e-16*(distance[0]*distance[1]*distance[2]
            +fabs(p[0]*distance[2])+fabs(p[1]*distance[0])+fabs(p[2]*distance[1]));
    }
    __device__ void solid(double h,double omega,double denominator) {
        height=h; angle=omega;
        double determinant=2.*area*h;
        angle_error=2.*fabs(determinant)*angle_error/(denominator*denominator+determinant*determinant);
    }
    __device__ void edge(int,double d,double a,double,double logarithm) {
        inside &= d>0.;
        absolute_terms+=fabs(d*logarithm);
        // Actual edge length is checked by the caller's geometry-independent
        // point/edge preflight below, not by a target/source index convention.
    }
    __device__ double error(double single) {
        if(fabs(height)<=64.*2.2204460492503131e-16*diameter && inside) flags|=1;
        if(!isfinite(single) || single<=0.) flags|=2;
        return (32.*2.2204460492503131e-16*(absolute_terms+fabs(height*angle))
                 +fabs(height)*angle_error)/(4.*3.14159265358979323846);
    }
};

extern "C" __global__ void paired_source_potential(
    const double* triangles,const double* normals,const double* outward,
    const double* tangents,const double* lengths,const double* areas,
    const double* centers,const double* signs,const double* diameters,
    const double* points,int targets,int panels,double image,double atol,double rtol,
    double* values,unsigned char* flags) {
    long long index=(long long)blockIdx.x*blockDim.x+threadIdx.x;
    if(index>=(long long)targets*panels) return;
    int target=index/panels, source=index%panels;
    double x[3],single[2],errors[2];
    unsigned char status=0;
    for(int side=0;side<2;++side) {
        for(int k=0;k<3;++k) x[k]=points[3*target+k];
        if(side) x[2]=2.*image-x[2];
        // Match the public CPU vertex/edge/near-plane rejection before using S.
        double h=0.;
        for(int k=0;k<3;++k) h+=(x[k]-triangles[9*source+k])*normals[3*source+k];
        for(int e=0;e<3;++e) {
            double d=0.,a=0.;
            for(int k=0;k<3;++k) {
                double r=x[k]-triangles[9*source+3*e+k];
                d-=r*outward[9*source+3*e+k];
                a+=r*tangents[9*source+3*e+k];
            }
            double l=lengths[3*source+e];
            double tolerance=64.*2.2204460492503131e-16;
            if(fabs(h)<=tolerance*diameters[source] && d*d+h*h<=(tolerance*l)*(tolerance*l)
                && a>=0. && a<=l) status|=1;
        }
        SourceGuard guard(diameters[source],areas[source]);
        double dbl,moment[3];
        rankine_sample_at_point_observed(triangles,normals,outward,tangents,lengths,
            areas,centers,signs,x,source,false,single[side],dbl,moment,guard);
        errors[side]=guard.error(single[side]);
        status|=guard.flags;
        if(!isfinite(errors[side]) || errors[side]>atol*diameters[source]+rtol*fabs(single[side])) status|=2;
    }
    double value=single[0]-single[1];
    double error=errors[0]+errors[1]+4.*2.2204460492503131e-16*(fabs(single[0])+fabs(single[1])+fabs(value));
    if(!isfinite(value) || !isfinite(error) || error>atol*diameters[source]+rtol*fabs(value)) status|=4;
    values[index]=value; flags[index]=status;
}
