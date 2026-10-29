// Copyright 2026 Alex Osborne. SPDX-License-Identifier: MIT
// Double-precision radix-2 convolution used by Bluestein's arbitrary-size DFT.
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL FP_CONTRACT OFF

double2 cmul(double2 a,double2 b) {
    return (double2)(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}
uint reversed(uint x,int bits) {
    x=((x&0x55555555U)<<1)|((x>>1)&0x55555555U);
    x=((x&0x33333333U)<<2)|((x>>2)&0x33333333U);
    x=((x&0x0f0f0f0fU)<<4)|((x>>4)&0x0f0f0f0fU);
    x=((x&0x00ff00ffU)<<8)|((x>>8)&0x00ff00ffU);
    x=(x<<16)|(x>>16);
    return x>>(32-bits);
}
__kernel void fft_stage(__global double2 *data,int m,int rows,int span,int inverse,
                        __global const double2 *twiddle) {
    size_t id=get_global_id(0);
    if(id>=(size_t)rows*m/2)return;
    int row=id/(m/2),b=id%(m/2),half_span=span/2,j=b%half_span;
    size_t a=(size_t)row*m+(b/half_span)*span+j,c=a+half_span;
    double2 t=twiddle[j*(m/span)];
    if(inverse)t.y=-t.y;
    double2 x=data[a],y=cmul(t,data[c]);
    data[a]=x+y;data[c]=x-y;
}
__kernel void start_real(__global const double *input,int n,int m,int rows,int bits,
                          __global const double2 *chirp,__global double2 *work) {
    size_t id=get_global_id(0);
    if(id>=(size_t)rows*m)return;
    int row=id/m,j=id%m;
    work[(size_t)row*m+reversed(j,bits)]=j<n ? input[(size_t)row*n+j]*chirp[j] : (double2)(0.,0.);
}
__kernel void start_inverse(__global const double2 *input,int n,int m,int rows,int bits,
                             __global const double2 *chirp,__global double2 *work) {
    size_t id=get_global_id(0);
    if(id>=(size_t)rows*m)return;
    int row=id/m,j=id%m,nr=n/2+1;
    double2 value=(double2)(0.,0.);
    if(j<n) {
        value=input[(size_t)row*nr+(j<=n/2 ? j : n-j)];
        if(j>n/2)value.y=-value.y;
        if(j==0 || (n%2==0 && j==n/2))value.y=0.;
        value.y=-value.y;
        value=cmul(value,chirp[j]);
    }
    work[(size_t)row*m+reversed(j,bits)]=value;
}
__kernel void convolve(__global const double2 *input,__global const double2 *convolution_values,
                        int m,int rows,int bits,__global double2 *output) {
    size_t id=get_global_id(0);
    if(id>=(size_t)rows*m)return;
    int row=id/m,j=id%m;
    output[(size_t)row*m+reversed(j,bits)]=cmul(input[id],convolution_values[j]);
}
__kernel void finish_forward(__global const double2 *work,__global const double2 *chirp,
                              int n,int m,int rows,__global double2 *output) {
    size_t id=get_global_id(0);
    int nr=n/2+1;
    if(id>=(size_t)rows*nr)return;
    int row=id/nr,j=id%nr;
    output[id]=cmul(chirp[j],work[(size_t)row*m+j])/(double)m;
}
__kernel void finish_selected(__global const double2 *work,__global const double2 *chirp,
                               int n,int m,int rows,int zero,__global const int *selected,
                               int count,__global double *output,__global double *initial,
                               __global double *premax) {
    size_t id=get_global_id(0);
    if(id>=(size_t)rows*count)return;
    int row=id/count,j=id%count,k=zero+selected[j];
    output[(size_t)j*rows+row]=cmul(chirp[k],work[(size_t)row*m+k]).x/((double)m*n);
    if(j==0) {
        initial[row]=cmul(chirp[zero],work[(size_t)row*m+zero]).x/((double)m*n);
        double peak=0.;
        for(int k0=0;k0<zero;++k0)
            peak=fmax(peak,fabs(cmul(chirp[k0],work[(size_t)row*m+k0]).x/((double)m*n)));
        premax[row]=peak;
    }
}
__kernel void multiply_spectrum(__global const double2 *H,int np,
        __global const double2 *velocity,int nr,int leading_zero,double scale,
        __global double2 *output) {
    size_t id=get_global_id(0);
    if(id>=(size_t)np*nr)return;
    int point=id/nr,f=id%nr;
    double2 sum=(double2)(0.,0.);
    if(!(leading_zero && f==0)) {
        for(int mode=0;mode<6;++mode) {
            double2 response=H[((size_t)(f-leading_zero)*np+point)*6+mode]*scale;
            response.y=-response.y;
            sum+=cmul(response,velocity[(size_t)mode*nr+f]);
        }
    }
    output[id]=sum;
}
