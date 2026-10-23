// Uses the SAME analytic sample primitive as energy and shape derivatives.
// This entry point has explicit nonself ownership for arbitrary graph targets.
template<bool gradients> struct SourceGuard {
    double diameter, area, height, angle, angle_error, absolute_terms;
    const double *edge_lengths,*outward,*normal;
    double gradient_absolute[3];
    bool inside;
    unsigned char flags;
    SourceGuard() = delete;
    __device__ SourceGuard(double d,double a,const double* l,const double* o,const double* n):
        diameter(d),area(a),height(0.),angle(0.),angle_error(0.),absolute_terms(0.),
        edge_lengths(l),outward(o),normal(n),inside(true),flags(0){
        for(int k=0;k<3;++k) gradient_absolute[k]=0.;
    }
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
    __device__ void edge(int edge,double d,double a,double,double logarithm) {
        inside &= d>0.;
        absolute_terms+=fabs(d*logarithm);
        const double tolerance=64.*2.2204460492503131e-16;
        const double l=edge_lengths[edge];
        if(fabs(height)<=tolerance*diameter && d*d+height*height<=(tolerance*l)*(tolerance*l)
            && a>=0. && a<=l) flags|=1;
        if(gradients) for(int k=0;k<3;++k) gradient_absolute[k]+=fabs(outward[3*edge+k]*logarithm);
    }
    __device__ double error(double single) {
        if(fabs(height)<=64.*2.2204460492503131e-16*diameter && inside) flags|=1;
        if(!isfinite(single) || single<=0.) flags|=2;
        return (32.*2.2204460492503131e-16*(absolute_terms+fabs(height*angle))
                 +fabs(height)*angle_error)/(4.*3.14159265358979323846);
    }
    __device__ double gradient_error(int k) {
        if(fabs(height)<=64.*2.2204460492503131e-16*diameter && inside) flags|=1;
        return (32.*2.2204460492503131e-16*(gradient_absolute[k]+fabs(angle*normal[k]))
                 +angle_error*fabs(normal[k]))/(4.*3.14159265358979323846);
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
        SourceGuard<false> guard(diameters[source],areas[source],lengths+3*source,
                                 outward+9*source,normals+3*source);
        double dbl,moment[3],gradient[3];
        rankine_sample_at_point_observed(triangles,normals,outward,tangents,lengths,
            areas,centers,signs,x,source,false,single[side],dbl,moment,gradient,guard);
        errors[side]=guard.error(single[side]);
        status|=guard.flags;
        if(!isfinite(errors[side]) || errors[side]>atol*diameters[source]+rtol*fabs(single[side])) status|=2;
    }
    double value=single[0]-single[1];
    double error=errors[0]+errors[1]+4.*2.2204460492503131e-16*(fabs(single[0])+fabs(single[1])+fabs(value));
    if(!isfinite(value) || !isfinite(error) || error>atol*diameters[source]+rtol*fabs(value)) status|=4;
    values[index]=value; flags[index]=status;
}

extern "C" __global__ void paired_source_gradient(
    const double* triangles,const double* normals,const double* outward,
    const double* tangents,const double* lengths,const double* areas,
    const double* centers,const double* signs,const double* diameters,
    const double* points,int targets,int panels,double image,double atol,double rtol,
    double* values,unsigned char* flags) {
    long long index=(long long)blockIdx.x*blockDim.x+threadIdx.x;
    if(index>=(long long)targets*panels) return;
    int target=index/panels,source=index%panels;
    double x[3],gradient[2][3],errors[2][3];
    unsigned char status=0;
    for(int side=0;side<2;++side) {
        for(int k=0;k<3;++k) x[k]=points[3*target+k];
        if(side) x[2]=2.*image-x[2];
        SourceGuard<true> guard(diameters[source],areas[source],lengths+3*source,
                                outward+9*source,normals+3*source);
        double single,dbl,moment[3];
        rankine_sample_at_point_observed(triangles,normals,outward,tangents,lengths,
            areas,centers,signs,x,source,false,single,dbl,moment,gradient[side],guard);
        for(int k=0;k<3;++k) {
            errors[side][k]=guard.gradient_error(k);
            if(!isfinite(gradient[side][k]) || !isfinite(errors[side][k])
                || errors[side][k]>atol+rtol*fabs(gradient[side][k])) status|=2;
        }
        status|=guard.flags;
    }
    for(int k=0;k<3;++k) {
        double sign=k==2?-1.:1.;
        double value=gradient[0][k]-sign*gradient[1][k];
        double error=errors[0][k]+errors[1][k]
            +4.*2.2204460492503131e-16*(fabs(gradient[0][k])+fabs(gradient[1][k])+fabs(value));
        if(!isfinite(value) || !isfinite(error) || error>atol+rtol*fabs(value)) status|=4;
        values[3*index+k]=value;
    }
    flags[index]=status;
}
