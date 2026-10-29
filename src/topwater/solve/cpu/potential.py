"""Potential-only host kernels for the shared reconstruction pipeline.

Same restricted domain, lookup tables, and potential-only formulas as the OpenCL
prototype. Frequency-independent Rankine potentials are prepared once per
target chunk, without gradients. See ../gpu/THIRD_PARTY.md and ../gpu/LICENSE.
"""
from contextlib import contextmanager
import numpy as np
from numba import njit, prange, get_num_threads, set_num_threads
from topwater.solve.capytaine import potential_inputs
from topwater.solve.transfer import Transfer


@njit(cache=True,inline='always')
def _lagrange(a,b,c,x):
    return ((x-a)*(x-b))/((c-a)*(c-b))


@njit(cache=True,inline='always')
def _wave(r,z,rr,zz,table,pi):
    nr,nz = len(rr),len(zz)
    if abs(z)<abs(zz[-1]) and r<rr[-1]:
        one = int(np.floor(nr*81/676.+.5))
        if r<1e-10:
            i=1
        elif r<1.:
            i=int(np.floor((np.log10(r)/10.+1.)*(one-1)+1.+.5))
        else:
            i=int(np.floor((r-1.)*(nr-one)/(rr[-1]-1.)+one+.5))
        dz=(np.log10(abs(zz[-1]))+10.)/nz
        jraw=(np.log10(max(abs(z),1e-300))+10.)/dz
        j=int(np.copysign(np.floor(abs(jraw)+.5),jraw))
        i=max(2,min(nr-1,i))-1
        j=max(2,min(nz-1,j))-1
        xl=(_lagrange(rr[i],rr[i+1],rr[i-1],r),
                     _lagrange(rr[i+1],rr[i-1],rr[i],r),
                     _lagrange(rr[i-1],rr[i],rr[i+1],r))
        zl=(_lagrange(zz[j],zz[j+1],zz[j-1],z),
                     _lagrange(zz[j+1],zz[j-1],zz[j],z),
                     _lagrange(zz[j-1],zz[j],zz[j+1],z))
        result=0j
        for a in range(3):
            row=0j
            for b in range(3):
                row+=table[i+a-1,j+b-1]*zl[b]
            result+=xl[a]*row
        return result
    invr=1./max(r,1e-10)
    invr1=1./np.hypot(r,z)
    amplitude=np.exp(z)*np.sqrt(2.*pi*invr)
    return 2.*complex(-amplitude*np.sin(r-pi/4.)+z*invr1**3-invr1,
                      amplitude*np.cos(r-pi/4.))


@njit(cache=True,parallel=True)
def _rankine(points,nodes,centers,normals,areas,radii,pi,result):
    for index in prange(len(points)*len(centers)):
        p,j=index//len(centers),index%len(centers)
        x=points[p]
        d=np.linalg.norm(x-centers[j])
        if d>7*radii[j]:
            value=areas[j]/d if d>1e-7 else 0.
        else:
            gz=np.dot(x-centers[j],normals[j])
            rr=np.empty(4)
            for a in range(4):rr[a]=np.linalg.norm(x-nodes[j,a])
            value=0.
            for a in range(4):
                b=(a+1)%4
                edge=nodes[j,b]-nodes[j,a]
                dk=np.linalg.norm(edge)
                if dk<float(np.float32(1e-3))*radii[j]:continue
                gy=np.dot(x-nodes[j,a],np.cross(normals[j],edge/dk))
                total=rr[a]+rr[b]
                ant=2*gy*dk
                dnt=total**2-dk**2+2*abs(gz)*total
                at=np.arctan(ant/dnt) if abs(gz)>=float(np.float32(1e-4))*radii[j] else 0.
                if abs(gy)<float(np.float32(1e-5)):
                    value-=2*at*abs(gz)
                else:
                    value+=gy*np.log((total+dk)/(total-dk))-2*at*abs(gz)
        result[p,j]=-value/(2*pi)
    return result


@njit(cache=True,inline='always')
def _potential_value(x,y,k,panel,centers,areas,quadrature,weights,rr,zz,table,pi,gamma,rankine):
    r=k*np.hypot(x-centers[panel,0],y-centers[panel,1])
    z=k*centers[panel,2]
    if abs(r)<float(np.float32(1e-10)) and abs(z)<float(np.float32(1e-10)):
        wave=k*areas[panel]*complex(1.-np.log(k*k*areas[panel]/pi)
                  +2.*(gamma-np.log(2.)),-2.*pi)
    else:
        wave=0j
        for q in range(weights.shape[1]):
            source=quadrature[panel,q]
            wave+=k*weights[panel,q]*_wave(k*np.hypot(x-source[0],y-source[1]),
                                          k*source[2],rr,zz,table,pi)
    return rankine-wave/(4.*pi)


@njit(cache=True,parallel=True)
def _potential(points,ks,centers,areas,quadrature,weights,rr,zz,table,pi,gamma,rankine,output):
    npnt,nf=len(points),len(centers)
    for index in prange(len(ks)*npnt*nf):
        panel=index%nf
        point=(index//nf)%npnt
        iw=index//(npnt*nf)
        output[iw,point,panel]=_potential_value(points[point,0],points[point,1],ks[iw],panel,
            centers,areas,quadrature,weights,rr,zz,table,pi,gamma,rankine[point,panel])


@njit(cache=True,parallel=True)
def _apply_sources(points,ks,centers,areas,quadrature,weights,rr,zz,table,pi,gamma,rankine,sources,output):
    for index in prange(len(ks)*len(points)):
        iw,point=index//len(points),index%len(points)
        for mode in range(sources.shape[2]):
            output[iw,point,mode]=0j
        for panel in range(len(centers)):
            value=_potential_value(points[point,0],points[point,1],ks[iw],panel,
                centers,areas,quadrature,weights,rr,zz,table,pi,gamma,rankine[point,panel])
            for mode in range(sources.shape[2]):
                output[iw,point,mode]+=value*sources[iw,panel,mode]


@njit(cache=True,parallel=True)
def _incident(xy,ks,directions,coefficients,output):
    for index in prange(len(ks)*len(xy)):
        iw,point=index//len(xy),index%len(xy)
        value=0j
        for carrier in range(len(directions)):
            plus,minus=coefficients[iw,carrier,0],coefficients[iw,carrier,1]
            if plus==0j and minus==0j:
                continue
            phase=ks[iw]*(xy[point,0]*directions[carrier,0]+xy[point,1]*directions[carrier,1])
            cosine,sine=np.cos(phase),np.sin(phase)
            value+=plus*complex(cosine,sine)
            value+=minus*complex(cosine,-sine)
        output[iw,point]=value


class CPUPotential:
    allocation_error = MemoryError

    def __init__(self, green_function, source_mesh, *, threads=None):
        self.inputs = potential_inputs(green_function, source_mesh)
        self.threads = min(8, get_num_threads()) if threads is None else threads
        self.table = np.ascontiguousarray(self.inputs.table[:, :, 0]+1j*self.inputs.table[:, :, 1])
        self.nf = len(self.inputs.centers)
        self._points = self._source_values = None
        self.output = self.rankine = self._rankine_storage = None

    @contextmanager
    def threaded(self):
        previous = get_num_threads()
        set_num_threads(self.threads)
        try:
            yield
        finally:
            set_num_threads(previous)

    def _output(self, shape):
        size = int(np.prod(shape))
        if self.output is None or self.output.size < size:
            self.output = np.empty(size, complex)
        return self.output[:size].reshape(shape)

    def _rankine_output(self, rows):
        size = rows*self.nf
        if self._rankine_storage is None or self._rankine_storage.size < size:
            self._rankine_storage = np.empty(size)
        return self._rankine_storage[:size].reshape(rows, self.nf)

    def reserve_radiation(self, rows, frequencies):
        self._output((frequencies, rows, 6))
        self._rankine_output(rows)

    @property
    def nbytes(self):
        return sum(value.nbytes for value in (self.output, self._rankine_storage) if value is not None)

    def reset_batch_buffers(self):
        self.output = self.rankine = self._rankine_storage = None
        self._points = self._source_values = None

    close = reset_batch_buffers

    def _arguments(self, points, wavenumbers):
        points = np.ascontiguousarray(points, dtype=float)
        ks = np.ascontiguousarray(wavenumbers, dtype=float)
        if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
                or not np.isfinite(points).all() or np.any(points[:, 2] != 0)):
            raise ValueError('Potential fields require finite nonempty (point,3) z=0 targets')
        if ks.ndim != 1 or not len(ks) or not np.isfinite(ks).all() or np.any(ks <= 0):
            raise ValueError('Potential fields require finite positive wavenumbers')
        return points, ks

    def _set_points(self, points):
        if self._points is not None and np.array_equal(points, self._points):
            return
        self._points = None
        g = self.inputs
        self.rankine = self._rankine_output(len(points))
        _rankine(points, g.nodes, g.centers, g.normals, g.areas, g.radii, g.pi, self.rankine)
        self._points = points.copy()

    def _kernel(self, function, points, ks, *tail):
        g = self.inputs
        function(points, ks, g.centers, g.areas, g.quadrature, g.weights,
                 g.rr, g.zz, self.table, g.pi, g.gamma, self.rankine, *tail)

    def evaluate(self, points, wavenumbers):
        points, ks = self._arguments(points, wavenumbers)
        output = np.empty((len(ks), len(points), self.nf), complex)
        with self.threaded():
            self._set_points(points)
            self._kernel(_potential, points, ks, output)
        return output

    def apply_sources(self, points, wavenumbers, sources):
        points, ks = self._arguments(points, wavenumbers)
        source = np.ascontiguousarray(sources, dtype=complex)
        if (source.ndim != 3 or source.shape[:2] != (len(ks), self.nf) or source.shape[2] < 1):
            raise ValueError('Sources must have (frequency,panel,mode) shape')
        if source.flags.writeable or source is not self._source_values:
            if not np.isfinite(source).all():
                raise ValueError('Sources must be finite')
            self._source_values = source
        output = self._output((len(ks), len(points), source.shape[2]))
        with self.threaded():
            self._set_points(points)
            self._kernel(_apply_sources, points, ks, source, output)
        return output

    def apply_sources_transfer(self, points, wavenumbers, sources):
        output = self.apply_sources(points, wavenumbers, sources)
        return Transfer(output, output.shape)

    def incident_spectrum(self, train, xy):
        from topwater.solve.validation import _finite
        points = _finite(xy, (None, 2), 'incident targets')
        coefficients = _finite(train.C, (None, 2, len(train.k)), 'incident coefficients', complex_value=True)
        theta = _finite(train.theta, (len(coefficients),), 'incident directions')
        ks = _finite(train.k, (len(train.k),), 'incident wavenumbers')
        directions = np.c_[np.cos(theta), np.sin(theta)]
        packed = np.ascontiguousarray(coefficients.transpose(2, 0, 1))
        output = self._output((len(ks), len(points)))
        with self.threaded():
            _incident(points, ks, directions, packed, output)
        return output
