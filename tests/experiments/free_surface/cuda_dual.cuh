// Shared forward directional arithmetic for the same analytic pair formula.
struct Dual {
    double v,d;
    __device__ Dual():v(0.),d(0.){}
    __device__ Dual(double a):v(a),d(0.){}
    __device__ Dual(double a,double b):v(a),d(b){}
    __device__ Dual& operator+=(Dual b){v+=b.v;d+=b.d;return *this;} __device__ Dual& operator-=(Dual b){v-=b.v;d-=b.d;return *this;}
};
__device__ Dual operator+(Dual a,Dual b){return Dual(a.v+b.v,a.d+b.d);}
__device__ Dual operator-(Dual a,Dual b){return Dual(a.v-b.v,a.d-b.d);}
__device__ Dual operator-(Dual a){return Dual(-a.v,-a.d);}
__device__ Dual operator*(Dual a,Dual b){return Dual(a.v*b.v,a.d*b.v+a.v*b.d);}
__device__ Dual operator/(Dual a,Dual b){return Dual(a.v/b.v,(a.d-(a.v/b.v)*b.d)/b.v);}
__device__ Dual sqrt(Dual a){double v=::sqrt(a.v);return Dual(v,a.d/(2.*v));}
__device__ Dual log1p(Dual a){return Dual(::log1p(a.v),a.d/(1.+a.v));}
__device__ Dual fabs(Dual a){return a.v<0.?-a:a;}
__device__ Dual fmax(Dual a,Dual b){return a.v>b.v?a:b;}
__device__ Dual atan2(Dual y,Dual x){return Dual(::atan2(y.v,x.v),(x.v*y.d-y.v*x.d)/(x.v*x.v+y.v*y.v));}
__device__ void atomicAdd(Dual* p,Dual a){atomicAdd(&p->v,a.v);atomicAdd(&p->d,a.d);}

