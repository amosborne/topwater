// Specific sample-local pullback, not a general AD tape. The shared owner
// evaluates every sqrt/log/angle exactly once and records only this sample.
struct SampleValues {
    double x[3],r[3][3],R[3],h,omega,solid_den,d[3],a[3],den[3],logarithm[3];
    __device__ void point(const double* p,const double v[3][3],const double* distance){
        for(int k=0;k<3;++k){x[k]=p[k];R[k]=distance[k];
            for(int t=0;t<3;++t)r[k][t]=v[k][t];}
    }
    __device__ void solid(double height,double angle,double denominator){
        h=height;omega=angle;solid_den=denominator;
    }
    __device__ void edge(int e,double distance,double along,double denominator,double logarithm_value){
        d[e]=distance;a[e]=along;den[e]=denominator;logarithm[e]=logarithm_value;
    }
};

template<int STRIDE>
__device__ __forceinline__ void add_vec(double out[9],const double* jac,int face,int offset,const double bar[3]){
    #pragma unroll
    for(int s=0;s<9;++s){
        #pragma unroll
        for(int k=0;k<3;++k)out[s]+=jac[(face*9+s)*STRIDE+offset+k]*bar[k];
    }
}
template<int STRIDE>
__device__ __forceinline__ void add_scalar(double out[9],const double* jac,int face,int offset,double bar){
    #pragma unroll
    for(int s=0;s<9;++s)out[s]+=jac[(face*9+s)*STRIDE+offset]*bar;
}

extern "C" __global__ void reverse_vertex_gradient(
    const double *tri0,const double *triJ,const double *normal0,const double *normalJ,
    const double *outward0,const double *outwardJ,const double *tangent0,const double *tangentJ,
    const double *length0,const double *lengthJ,const double *area0,const double *areaJ,
    const double *center0,const double *centerJ,const double *sign0,const double *signJ,
    const double *gradient0,const double *gradientJ,const double *points0,const double *pointsJ,
    const double *weights,const long long *faces,const double *curlphi0,const double *curlphiJ,
    const double *g0,const double *gJ,const double *phi,const double *aux,
    int nf,int nv,int nq,double density,double *out
#ifdef PAIR_DIAGNOSTIC
    ,const long long *selected,int selected_count
#endif
    ){
    int work=blockDim.x*blockIdx.x+threadIdx.x;
#ifdef PAIR_DIAGNOSTIC
    if(work>=selected_count)return;
    int slot=work,pair=selected[slot];
    double pair_value=0.;
#else
    int pair=work;
    if(pair>=nf*nf)return;
#endif
    int i=pair/nf,j=pair%nf;
    double zi=aux[i],sign=sign0[j],g[3],A=-zi*aux[j],sum_phi=0.;
    #pragma unroll
    for(int k=0;k<3;++k){g[k]=g0[j*3+k];A+=curlphi0[i*3+k]*curlphi0[j*3+k];}
    A*=.5;
    for(int c=0;c<3;++c)sum_phi+=phi[faces[j*3+c]];
    double target[9]={0.},source[9]={0.},barA=0.,barg[3]={0.};
    double target_area_bar=0.;
    for(int q=0;q<nq;++q){
        SampleValues p;double single,dbl,moment[3];
        energy_sample_observed<double,const double*,SampleValues>(tri0,normal0,
            outward0,tangent0,length0,area0,center0,sign0,points0,
            i,j,nq,q,single,dbl,moment,p);
        double mixed=dbl*(sum_phi/3.);
        for(int k=0;k<3;++k)mixed+=g[k]*moment[k];
        double weight=density*area0[i]*weights[q];
        double kappa=weight/(4.*3.1415926535897932384626433832795);
        target_area_bar+=density*weights[q]*(A*single+zi*mixed);
#ifdef PAIR_DIAGNOSTIC
        pair_value+=weight*(A*single+zi*mixed);
#endif
        barA+=weight*single;
        for(int k=0;k<3;++k)barg[k]+=weight*zi*moment[k];
        double B=sum_phi/3.,ng=0.,sum_og_log=0.;
        for(int k=0;k<3;++k){
            B+=g[k]*(p.x[k]-center0[j*3+k]-p.h*normal0[j*3+k]);
            ng+=normal0[j*3+k]*g[k];
        }
        for(int e=0;e<3;++e){
            double og=0.;for(int k=0;k<3;++k)og+=outward0[j*9+e*3+k]*g[k];
            sum_og_log+=og*p.logarithm[e];
        }
        double bar_h=kappa*(-zi*sign*sum_og_log-A*p.omega-zi*sign*p.omega*ng);
        double bar_omega=kappa*(zi*sign*B-A*p.h);
        double bar_r[3][3]={{0.}},bar_R[3]={0.},bar_x[3],bar_n[3],bar_center[3];
        for(int k=0;k<3;++k){
            bar_x[k]=kappa*zi*sign*p.omega*g[k];
            bar_center[k]=-bar_x[k];bar_n[k]=-p.h*bar_x[k];
        }
        // Edge pullbacks are reduced into source coordinate bars immediately.
        for(int e=0;e<3;++e){
            int next=(e+1)%3;double og=0.;
            for(int k=0;k<3;++k)og+=outward0[j*9+e*3+k]*g[k];
            double l=length0[j*3+e],a=p.a[e],end=l-a,d=p.d[e];
            double perp=d*d+p.h*p.h,D=p.den[e],u=2.*l/D;
            double bar_log=kappa*(A*d-zi*p.h*sign*og);
            double bar_u=bar_log/(1.+u),bar_l=bar_u*2./D,bar_D=-bar_u*u/D;
            double s0=p.R[e]+fabs(a),s1=p.R[next]+fabs(end);
            double bar_p=bar_D*(1./s0+1./s1);
            double b0=-bar_D*perp/(s0*s0),b1=-bar_D*perp/(s1*s1);
            bar_R[e]+=b0;bar_R[next]+=b1;
            // Match Dual fabs(0) and strict fmax tie branches.
            double abs_a=a<0.?-1.:1.,abs_end=end<0.?-1.:1.;
            double bar_a=b0*abs_a-b1*abs_end;
            bar_l+=b1*abs_end;
            if(-a>0.)bar_a-=2.*bar_D;
            if(a-l>0.){bar_a+=2.*bar_D;bar_l-=2.*bar_D;}
            double bar_d=kappa*A*p.logarithm[e]+2.*d*bar_p;
            bar_h+=2.*p.h*bar_p;
            double bar_o[3],bar_t[3];
            for(int k=0;k<3;++k){
                bar_r[e][k]+=-bar_d*outward0[j*9+e*3+k]+bar_a*tangent0[j*9+e*3+k];
                bar_o[k]=-bar_d*p.r[e][k]-kappa*zi*p.h*sign*p.logarithm[e]*g[k];
                bar_t[k]=bar_a*p.r[e][k];
            }
            add_vec<9>(source,outwardJ,j,e*3,bar_o);
            add_vec<9>(source,tangentJ,j,e*3,bar_t);
            add_scalar<3>(source,lengthJ,j,e,bar_l);
        }
        if(i!=j){
            double N=2.*area0[j]*p.h,D=p.solid_den;
            double denom=D*D+N*N;
            double bar_N=bar_omega*2.*D/denom,bar_D=-bar_omega*2.*N/denom;
            add_scalar<1>(source,areaJ,j,0,bar_N*2.*p.h);
            bar_h+=bar_N*2.*area0[j];
            for(int k=0;k<3;++k)bar_R[k]+=bar_D*p.R[(k+1)%3]*p.R[(k+2)%3];
            for(int e=0;e<3;++e){
                int next=(e+1)%3,other=(e+2)%3;double dot=0.;
                for(int k=0;k<3;++k){
                    dot+=p.r[e][k]*p.r[next][k];
                    bar_r[e][k]+=bar_D*p.R[other]*p.r[next][k];
                    bar_r[next][k]+=bar_D*p.R[other]*p.r[e][k];
                }
                bar_R[other]+=bar_D*dot;
            }
            for(int k=0;k<3;++k){bar_r[0][k]+=bar_h*normal0[j*3+k];bar_n[k]+=bar_h*p.r[0][k];}
        }
        for(int e=0;e<3;++e){
            double bar_tri[3];
            for(int k=0;k<3;++k){
                bar_r[e][k]+=bar_R[e]*p.r[e][k]/p.R[e];
                bar_x[k]+=bar_r[e][k];bar_tri[k]=-bar_r[e][k];
            }
            add_vec<9>(source,triJ,j,e*3,bar_tri);
        }
        add_vec<3>(source,normalJ,j,0,bar_n);
        add_vec<3>(source,centerJ,j,0,bar_center);
        // pointsJ stride depends on nq; keep the same owner/layout explicitly.
        #pragma unroll
        for(int s=0;s<9;++s){
            #pragma unroll
            for(int k=0;k<3;++k)target[s]+=pointsJ[((i*9+s)*nq+q)*3+k]*bar_x[k];
        }
    }
    add_scalar<1>(target,areaJ,i,0,target_area_bar+(i==j?density*zi*sum_phi/6.:0.));
    double bar_ci[3],bar_cj[3];
    for(int k=0;k<3;++k){bar_ci[k]=.5*barA*curlphi0[j*3+k];bar_cj[k]=.5*barA*curlphi0[i*3+k];}
    add_vec<3>(target,curlphiJ,i,0,bar_ci);
    add_vec<3>(source,curlphiJ,j,0,bar_cj);
    add_vec<3>(source,gJ,j,0,barg);
    #pragma unroll
    for(int s=0;s<9;++s){
#ifdef PAIR_DIAGNOSTIC
        out[slot*19+s]=target[s];out[slot*19+9+s]=source[s];
#else
        atomicAdd(out+faces[i*3+s/3]*3+s%3,target[s]);
        atomicAdd(out+faces[j*3+s/3]*3+s%3,source[s]);
#endif
    }
#ifdef PAIR_DIAGNOSTIC
    out[slot*19+18]=pair_value+(i==j?density*zi*area0[i]*sum_phi/6.:0.);
#endif
}
