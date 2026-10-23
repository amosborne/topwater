// Fixed-width forward derivatives share each primal elementary evaluation.
template<int N> struct WideDual {
    double v,d[N];
    __device__ WideDual():v(0.){for(int k=0;k<N;++k)d[k]=0.;}
    __device__ WideDual(double a):v(a){for(int k=0;k<N;++k)d[k]=0.;}
    __device__ WideDual& operator+=(WideDual b){v+=b.v;for(int k=0;k<N;++k)d[k]+=b.d[k];return *this;}
    __device__ WideDual& operator-=(WideDual b){v-=b.v;for(int k=0;k<N;++k)d[k]-=b.d[k];return *this;}
};
template<int N> __device__ WideDual<N> operator+(WideDual<N> a,WideDual<N> b){
    WideDual<N> r(a.v+b.v);for(int k=0;k<N;++k)r.d[k]=a.d[k]+b.d[k];return r;}
template<int N> __device__ WideDual<N> operator-(WideDual<N> a,WideDual<N> b){
    WideDual<N> r(a.v-b.v);for(int k=0;k<N;++k)r.d[k]=a.d[k]-b.d[k];return r;}
template<int N> __device__ WideDual<N> operator-(WideDual<N> a){
    WideDual<N> r(-a.v);for(int k=0;k<N;++k)r.d[k]=-a.d[k];return r;}
template<int N> __device__ WideDual<N> operator*(WideDual<N> a,WideDual<N> b){
    WideDual<N> r(a.v*b.v);for(int k=0;k<N;++k)r.d[k]=a.d[k]*b.v+a.v*b.d[k];return r;}
template<int N> __device__ WideDual<N> operator/(WideDual<N> a,WideDual<N> b){
    WideDual<N> r(a.v/b.v);for(int k=0;k<N;++k)r.d[k]=(a.d[k]-r.v*b.d[k])/b.v;return r;}
template<int N> __device__ WideDual<N> sqrt(WideDual<N> a){
    WideDual<N> r(::sqrt(a.v));for(int k=0;k<N;++k)r.d[k]=a.d[k]/(2.*r.v);return r;}
template<int N> __device__ WideDual<N> log1p(WideDual<N> a){
    WideDual<N> r(::log1p(a.v));for(int k=0;k<N;++k)r.d[k]=a.d[k]/(1.+a.v);return r;}
template<int N> __device__ WideDual<N> fabs(WideDual<N> a){return a.v<0.?-a:a;}
template<int N> __device__ WideDual<N> fmax(WideDual<N> a,WideDual<N> b){return a.v>b.v?a:b;}
template<int N> __device__ WideDual<N> atan2(WideDual<N> y,WideDual<N> x){
    WideDual<N> r(::atan2(y.v,x.v));double denominator=x.v*x.v+y.v*y.v;
    for(int k=0;k<N;++k)r.d[k]=(x.v*y.d[k]-y.v*x.d[k])/denominator;return r;}
// Mixed scalar operations cannot deduce N through an implicit conversion.
template<int N> __device__ WideDual<N> operator+(double a,WideDual<N> b){return WideDual<N>(a)+b;}
template<int N> __device__ WideDual<N> operator+(WideDual<N> a,double b){return a+WideDual<N>(b);}
template<int N> __device__ WideDual<N> operator-(double a,WideDual<N> b){return WideDual<N>(a)-b;}
template<int N> __device__ WideDual<N> operator-(WideDual<N> a,double b){return a-WideDual<N>(b);}
template<int N> __device__ WideDual<N> operator*(double a,WideDual<N> b){return WideDual<N>(a)*b;}
template<int N> __device__ WideDual<N> operator*(WideDual<N> a,double b){return a*WideDual<N>(b);}
template<int N> __device__ WideDual<N> operator/(double a,WideDual<N> b){return WideDual<N>(a)/b;}
template<int N> __device__ WideDual<N> operator/(WideDual<N> a,double b){return a/WideDual<N>(b);}
template<int N> __device__ WideDual<N> fmax(WideDual<N> a,double b){return fmax(a,WideDual<N>(b));}
template<int N> __device__ WideDual<N> fmax(double a,WideDual<N> b){return fmax(WideDual<N>(a),b);}
