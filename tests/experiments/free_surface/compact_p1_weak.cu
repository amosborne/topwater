// Whole-hat weak |D| kernel. Near/touching supports retain the complete
// gradient--single-layer form. Only geometrically separated WHOLE hats use
// the cancellation-free inverse-cube integral.
struct WeakSampleObserver {
    double absolute_terms;
    __device__ WeakSampleObserver():absolute_terms(0.){}
    __device__ void point(const double*,const double[3][3],const double*){}
    __device__ void solid(double,double,double){}
    __device__ void edge(int,double d,double,double,double logarithm){absolute_terms+=fabs(d*logarithm);}
};

__device__ double accumulation_error(double count) {
    const double product=count*2.2204460492503131e-16;
    return product<1.?product/(1.-product):nan("");
}

__device__ double near_hat(
    const double* tri,const double* normal,const double* outward,const double* tangent,
    const double* length,const double* area,const double* center,const double* sign,
    const double* gradient,const double* points,const double* weights,
    const int* cells,const int* orientation,int nq,int dx,int dy,double hx,double hy,double& error) {
    double total=0.,absolute=0.,count=0.;error=0.;
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) {
        const double dot=gradient[2*i]*gradient[2*j]+gradient[2*i+1]*gradient[2*j+1];
        if(dot==0.) continue;
        const bool same=cells[2*i]+dx==cells[2*j]
            && cells[2*i+1]+dy==cells[2*j+1] && orientation[i]==orientation[j];
        for(int q=0;q<nq;++q) {
            double x[3]={points[(i*nq+q)*3]+dx*hx,points[(i*nq+q)*3+1]+dy*hy,0.};
            double single,dbl,moment[3],gs[3];
            WeakSampleObserver observer;
            rankine_sample_at_point_observed(tri,normal,outward,tangent,length,area,center,
                sign,x,j,same,single,dbl,moment,gs,observer);
            if(!isfinite(single) || single<=0.) return nan("");
            const double weight=2.*dot*area[i]*weights[q];
            total+=weight*single;absolute+=fabs(weight*single);
            count+=1.;
            // All targets/sources are coplanar: height*angle is exactly zero,
            // including owned self. Retain inner S cancellation separately.
            error+=fabs(weight)*32.*2.2204460492503131e-16*observer.absolute_terms/(4.*3.14159265358979323846);
        }
    }
    error+=(64.*2.2204460492503131e-16+accumulation_error(count))*absolute;
    return total;
}

extern "C" __global__ void weak_hat_lags(
    const double* tri,const double* normal,const double* outward,const double* tangent,
    const double* length,const double* area,const double* center,const double* sign,
    const double* gradient,const double* points,const double* weights,
    const int* cells,const int* orientation,int nq,
    const double* hat_points,const double* hat_weights,int nh,
    const int* lags,int count,double hx,double hy,double* value,double* asymmetry,double* roundoff) {
    int row=blockIdx.x*blockDim.x+threadIdx.x;
    if(row>=count) return;
    int dx=lags[2*row],dy=lags[2*row+1];
    double sx=fmax(double(abs(dx)-2),0.)*hx,sy=fmax(double(abs(dy)-2),0.)*hy;
    const bool separated=sx*sx+sy*sy>=hx*hx+hy*hy;
    if(separated) {
        double sum=0.;
        for(int a=0;a<nh;++a) for(int b=0;b<nh;++b) {
            double x=dx*hx+hat_points[2*a]-hat_points[2*b];
            double y=dy*hy+hat_points[2*a+1]-hat_points[2*b+1];
            double r2=x*x+y*y;
            sum+=hat_weights[a]*hat_weights[b]/(r2*sqrt(r2));
        }
        value[row]=-sum/(2.*3.14159265358979323846);
        asymmetry[row]=0.;
        roundoff[row]=(64.*2.2204460492503131e-16+accumulation_error(double(nh)*nh))*fabs(value[row]);
    } else {
        double ea,eb;
        double a=near_hat(tri,normal,outward,tangent,length,area,center,sign,gradient,
            points,weights,cells,orientation,nq,dx,dy,hx,hy,ea);
        double b=near_hat(tri,normal,outward,tangent,length,area,center,sign,
            gradient,points,weights,cells,orientation,nq,-dx,-dy,hx,hy,eb);
        // Symmetric target/source quadrature, not an eigenvalue/row repair.
        value[row]=(a+b)*.5;
        asymmetry[row]=fabs(a-b);
        roundoff[row]=(ea+eb)*.5;
    }
}
