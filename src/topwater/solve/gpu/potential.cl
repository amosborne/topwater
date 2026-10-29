// Potential-only translation of Capytaine 3.0.0 libDelhommeau routines.
// Copyright 2026 Capytaine developers; modifications Copyright 2026 Alex Osborne.
// SPDX-License-Identifier: Apache-2.0
// Derived from Green_Rankine.f90, Green_wave.f90 and Delhommeau_integrals.f90.
// See THIRD_PARTY.md. Preserve the reference's quadrature and runtime constants.
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL FP_CONTRACT OFF

double rankine(double3 x, int panel, __global const double *nodes,
               __global const double *centers, __global const double *normals,
               __global const double *areas, __global const double *radii) {
    double3 center = vload3(panel, centers);
    double distance = length(x-center), radius = radii[panel];
    if (distance > 7.0*radius)
        return distance > 1e-7 ? areas[panel]/distance : 0.0;
    double3 normal = vload3(panel, normals);
    double gz = dot(x-center, normal), rr[4];
    for (int a=0; a<4; ++a) rr[a] = length(x-vload3(4*panel+a, nodes));
    double result = 0.0;
    for (int a=0; a<4; ++a) {
        int b = (a+1)%4;
        double3 start = vload3(4*panel+a, nodes);
        double3 edge = vload3(4*panel+b, nodes)-start;
        double dk = length(edge);
        if (dk < (double)(float)1e-3*radius) continue;
        double gy = dot(x-start, cross(normal, edge/dk));
        double sum = rr[a]+rr[b];
        double ant = 2.0*gy*dk;
        double dnt = sum*sum-dk*dk+2.0*fabs(gz)*sum;
        double at = fabs(gz) >= (double)(float)1e-4*radius ? atan(ant/dnt) : 0.0;
        if (fabs(gy) < (double)(float)1e-5)
            result -= 2.0*at*fabs(gz);
        else
            result += gy*log((sum+dk)/(sum-dk))-2.0*at*fabs(gz);
    }
    return result;
}

__kernel void rankine_matrix(__global const double *points, int np, int nf,
        __global const double *nodes, __global const double *centers,
        __global const double *normals, __global const double *areas,
        __global const double *radii, double pi, __global double *output) {
    int index = get_global_id(0);
    if (index >= np*nf) return;
    int point = index/nf, panel = index%nf;
    // At z=0 the direct and reflected Rankine potentials are equal.
    output[index] = -rankine(vload3(point, points), panel, nodes, centers,
                            normals, areas, radii)/(2.0*pi);
}

double lagrange(double a, double b, double c, double x) {
    return ((x-a)*(x-b))/((c-a)*(c-b));
}

double2 wave(double r, double z, int nr, int nz, double pi,
        __global const double *rr, __global const double *zz,
        __global const double2 *table) {
    if (fabs(z) < fabs(zz[nz-1]) && r < rr[nr-1]) {
        // Scaled Nemoh3 indexing, including clamping at the table boundaries.
        int one = (int)floor((double)(nr*81)/676.0+.5);
        int i;
        if (r < 1e-10) i = 1;
        else if (r < 1.0) i = (int)floor((log10(r)/10.0+1.0)*(one-1)+1.0+.5);
        else i = (int)floor((r-1.0)*(nr-one)/(rr[nr-1]-1.0)+one+.5);
        double dz = (log10(fabs(zz[nz-1]))+10.0)/nz;
        double jraw = (log10(fmax(fabs(z), 1e-300))+10.0)/dz;
        int j = (int)copysign(floor(fabs(jraw)+.5), jraw);
        i = clamp(i, 2, nr-1)-1; j = clamp(j, 2, nz-1)-1;
        double xl[3], zl[3];
        xl[0] = lagrange(rr[i], rr[i+1], rr[i-1], r);
        xl[1] = lagrange(rr[i+1], rr[i-1], rr[i], r);
        xl[2] = lagrange(rr[i-1], rr[i], rr[i+1], r);
        zl[0] = lagrange(zz[j], zz[j+1], zz[j-1], z);
        zl[1] = lagrange(zz[j+1], zz[j-1], zz[j], z);
        zl[2] = lagrange(zz[j-1], zz[j], zz[j+1], z);
        double2 result = (double2)(0.0, 0.0);
        for (int a=0; a<3; ++a) {
            double2 row = (double2)(0.0, 0.0);
            for (int b=0; b<3; ++b) row += table[(i+a-1)*nz+j+b-1]*zl[b];
            result += xl[a]*row;
        }
        return result;
    }
    double invr = 1.0/fmax(r, 1e-10), invr1 = 1.0/hypot(r,z);
    double amplitude = exp(z)*sqrt(2.0*pi*invr);
    return 2.0*(double2)(-amplitude*sin(r-pi/4.0)+z*invr1*invr1*invr1-invr1,
                        amplitude*cos(r-pi/4.0));
}

__kernel void potential_matrix(__global const double *points, int np, int nf,
        __global const double *centers, __global const double *areas,
        int nq, __global const double *quadrature, __global const double *weights,
        __global const double *ks, int nk, int nr, int nz,
        __global const double *rr, __global const double *zz,
        __global const double2 *table, double pi, double gamma,
        __global const double *rankine, __global double2 *output) {
    size_t index = get_global_id(0);
    if (index >= (size_t)nk*np*nf) return;
    int panel = index%nf, point = (index/nf)%np, iw = index/((size_t)np*nf);
    double k = ks[iw];
    double3 x = vload3(point, points), center = vload3(panel, centers);
    double r = k*length(x.xy-center.xy), z = k*center.z;
    double2 sum = (double2)(0.0, 0.0);
    if (fabs(r) < (double)(float)1e-10 && fabs(z) < (double)(float)1e-10) {
        sum = k*areas[panel]*(double2)(
            1.0-log(k*k*areas[panel]/pi)+2.0*(gamma-log(2.0)), -2.0*pi);
    } else {
        for (int q=0; q<nq; ++q) {
            double3 source = vload3(panel*nq+q, quadrature);
            sum += k*weights[panel*nq+q]*wave(
                k*length(x.xy-source.xy), k*source.z, nr,nz,pi,rr,zz,table);
        }
    }
    output[index] = (double2)(rankine[point*nf+panel],0.0)-sum/(4.0*pi);
}

double2 multiply_complex(double2 a,double2 b) {
    return (double2)(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}

__kernel void apply_sources(__global const double2 *potential,
        __global const double2 *sources,int nf,int np,int nm,int nk,int offset,
        int target_np,int target_point,int target_frequency,
        __global double2 *output) {
    size_t index=get_global_id(0);
    if(index>=(size_t)nk*np*nm)return;
    int mode=index%nm,point=(index/nm)%np,iw=index/((size_t)np*nm);
    double2 value=(double2)(0.0,0.0);
    for(int j=0;j<nf;++j)
        value+=multiply_complex(potential[((size_t)iw*np+point)*nf+j],
                                sources[((size_t)(iw+offset)*nf+j)*nm+mode]);
    output[((size_t)(iw+target_frequency)*target_np+point+target_point)*nm+mode]=value;
}

// Exact physical plane-wave sum; coefficients already contain quadrature weights.
__kernel void incident_spectrum(__global const double *xy,int np,
        __global const double *ks,int nk,__global const double *directions,int nc,
        __global const double2 *coefficients,__global double2 *output) {
    int index=get_global_id(0);
    if(index>=nk*np)return;
    int point=index%np,iw=index/np;
    double x=xy[2*point],y=xy[2*point+1];
    double2 value=(double2)(0.0,0.0);
    for(int carrier=0;carrier<nc;++carrier) {
        double2 plus=coefficients[(iw*nc+carrier)*2];
        double2 minus=coefficients[(iw*nc+carrier)*2+1];
        if(plus.x==0.0 && plus.y==0.0 && minus.x==0.0 && minus.y==0.0)continue;
        double phase=ks[iw]*(x*directions[2*carrier]+y*directions[2*carrier+1]);
        double cosine;
        double sine=sincos(phase,&cosine);
        value+=multiply_complex(plus,(double2)(cosine,sine));
        value+=multiply_complex(minus,(double2)(cosine,-sine));
    }
    output[index]=value;
}
