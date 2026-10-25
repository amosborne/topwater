"""One exact dispersive signed-kappa train shared by forcing, calm bounds and video."""
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import czt


def window(tau):
    r = np.abs(np.asarray(tau))
    q = np.clip((r-15)/10,0,1)
    return 1-10*q**3+15*q**4-6*q**5


def window_spectrum(nu, nodes=256):
    """G(nu)=(1/2pi) integral W(tau)exp(-i nu tau)dtau; no L1/L2 rescaling."""
    x,w = np.polynomial.legendre.leggauss(nodes)
    tau = 20+5*x
    weights = 5*w*window(tau)/np.pi
    flat = np.asarray(nu).reshape(-1)
    result = np.empty_like(flat)
    for start in range(0,len(flat),2048):
        v = flat[start:start+2048]
        result[start:start+2048] = 15/np.pi*np.sinc(15*v/np.pi)+np.cos(v[:,None]*tau)@weights
    return result.reshape(np.shape(nu))


class CalmBuoyTrain:
    def __init__(self, path):
        with np.load(path) as stored:
            self.data = {key:stored[key].copy() for key in stored.files}
        self.omega = self.data['omega_rad_s']
        self.k = self.data['k_abs_m_inv']
        self.C = self.data['C']
        self.theta = self.data['directions_body_rad']
        assert self.C.shape == (384,2,len(self.omega))

    @staticmethod
    def construct(carriers, case, controls, g, beam):
        a = carriers['amplitudes']
        center_time = float(controls['center_group_time_s'])
        assert center_time in [30.,40.] and controls['plateau_halfwidth_s'] == 15
        assert controls['support_halfwidth_s'] == 25 and controls['sideband_omega_min_rad_s'] == 0
        assert controls['sideband_omega_step_rad_s'] == .005 and controls['signed_kappa'] == [1,-1]
        w0 = carriers['angular_frequencies']
        phi = carriers['phases_body']
        cg = g/(2*w0)
        k0 = w0**2/g
        maximum = 2*np.pi*min(.6,np.sqrt(g/(2*np.pi*2*beam)))
        omega = np.arange(0.,maximum,.005)
        omega = np.r_[omega,maximum]
        d = np.diff(omega)
        weights = np.r_[d[0]/2,(d[:-1]+d[1:])/2,d[-1]/2]
        k = omega**2/g
        signs = np.array([1.,-1.])
        nu = cg[:,None,None]*(signs[None,:,None]*k[None,None,:]-k0[:,None,None])
        G = window_spectrum(nu,controls['transition_gauss_nodes'])
        jacobian = (2*cg[:,None]*omega[None,:]/g)*weights[None,:]
        C = a[:,None,None]*G*np.exp(1j*(phi[:,None,None]+center_time*nu))*jacobian[:,None,:]
        # Componentwise spectral norm per crest length: a^2*cg*integral |G|^2 dnu.
        x,w = np.polynomial.legendre.leggauss(controls['transition_gauss_nodes'])
        full_G2_integral = (30+2*5*np.dot(w,window(20+5*x)**2))/(2*np.pi)
        full_norm = a*a*cg*full_G2_integral
        norm_weights = a[:,None,None]**2*cg[:,None,None]*G*G*jacobian[:,None,:]
        retained_norm = norm_weights.sum(axis=(1,2))
        return dict(omega_rad_s=omega,k_abs_m_inv=k,C=C,signs=signs,
            directions_body_rad=carriers['directions_body'],
            directions_global_rad=carriers['directions'],component_id=carriers['component_id'],
            band_index=carriers['band_index'],carrier_amplitudes_m=a,carrier_omega_rad_s=w0,
            carrier_phases_body_rad=phi,carrier_phases_global_rad=carriers['phases'],
            global_from_boat_rotation=np.asarray(case['global_from_boat_rotation']),
            global_from_boat_translation_m=np.asarray(case['global_from_boat_translation_m']),
            full_component_spectral_norm=full_norm,retained_component_spectral_norm=retained_norm,
            spectral_norm_weights=norm_weights,
            represented_reverse_component_norm=a*a*cg*np.sum(G[:,1]**2*jacobian,axis=1),
            window_center_group_time_s=center_time,gravity_m_s2=g)

    @staticmethod
    def static_excitation(table):
        triangles = np.asarray(table['triangles'])
        cg = np.asarray(table['cg']).reshape(3)
        pressure = float(np.asarray(table['rho']).item())*float(np.asarray(table['g']).item())
        area_normal = np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])/2
        panel_force = -pressure*area_normal
        H0 = np.r_[panel_force.sum(axis=0),np.cross(triangles.mean(axis=1)-cg,panel_force).sum(axis=0)]
        Kheave = np.asarray(table['K'])[:,2]
        np.testing.assert_allclose(H0,Kheave,rtol=0,atol=1e-8)
        return H0,float(np.max(abs(H0-Kheave)))

    def force(self, table, times):
        omega = np.asarray(table['omega_rad_s']).reshape(-1)
        theta = np.asarray(table['heading_rad']).reshape(-1)
        H = np.asarray(table['H_complex'])
        assert H.shape == (34,25,6) and theta[0] == 0 and abs(theta[-1]-2*np.pi) < 1e-12
        assert np.array_equal(H[:,0],H[:,-1])
        assert self.omega.min() == 0 and omega.min() == .1 and self.omega.max() <= omega.max()
        interp = RegularGridInterpolator((omega,theta),H,method='linear',bounds_error=True)
        headings = np.mod(self.theta[:,None]+np.array([0.,np.pi])[None,:],2*np.pi)
        W = np.broadcast_to(self.omega,self.C.shape)
        TH = np.broadcast_to(headings[:,:,None],self.C.shape)
        # Only lowfrequency endpoint evaluation uses .1; physical C is NOT clipped.
        actualH = interp(np.stack([np.maximum(W,.1),TH],axis=-1))
        H0,_ = self.static_excitation(table)
        low = self.omega < .1
        actualH[:,:,low,:] = H0+(self.omega[low][None,None,:,None]/.1)**2*(actualH[:,:,low,:]-H0)
        assert np.all(self.C[:,:,0] == 0),'Jacobian atOmega0 must give no DC coefficient'
        spectral_force = np.sum(self.C[:,:,:,None]*actualH,axis=(0,1))
        force = (np.exp(-1j*np.asarray(times)[:,None]*self.omega)@spectral_force).real
        carrier_queries = np.c_[self.data['carrier_omega_rad_s'],np.mod(self.theta,2*np.pi)]
        carrierH = interp(carrier_queries)
        scale = np.sum(self.data['carrier_amplitudes_m'][:,None]*abs(carrierH),axis=0)
        return force,scale,spectral_force

    def elevation(self, xy, times):
        """xy is NATIVE mean-boat XY; global points require R.T@(global-t) once."""
        xy = np.atleast_2d(xy)
        times = np.atleast_1d(times)
        out = np.zeros((len(times),len(xy)))
        temporal = np.exp(-1j*times[:,None]*self.omega)
        for j in np.flatnonzero(self.data['carrier_amplitudes_m']):
            s = xy@np.array([np.cos(self.theta[j]),np.sin(self.theta[j])])
            spatial = (self.C[j,0,:,None]*np.exp(1j*self.k[:,None]*s)+
                       self.C[j,1,:,None]*np.exp(-1j*self.k[:,None]*s))
            out += (temporal@spatial).real
        return out

    def profiles(self, component_ids, s_grid, times):
        """Video cache API: float64[component,time,s], exact spectral sum via CZT.

        s is NATIVE mean-boat projection along the carrier's body direction.
        Global points require R.T@(global-t) once, not just changing theta.
        Call in small component batches; linearly interpolate only the cached
        spatial profiles for video. Never use cache interpolation for forcing.
        """
        ids = np.asarray(component_ids,dtype=int)
        s = np.asarray(s_grid)
        t = np.asarray(times)
        assert len(t) >= 2 and np.max(abs(np.diff(t)-np.diff(t)[0])) < 1e-10
        assert np.max(abs(np.diff(self.omega[:-1])-.005)) < 1e-12
        B = (self.C[ids,0,:,None]*np.exp(1j*self.k[:,None]*s)+
             self.C[ids,1,:,None]*np.exp(-1j*self.k[:,None]*s))
        values = czt(B[:,:-1,:],m=len(t),w=np.exp(-1j*.005*(t[1]-t[0])),
                     a=np.exp(1j*.005*t[0]),axis=1)
        values *= np.exp(-1j*self.omega[0]*t)[None,:,None]
        values += B[:,-1:, :]*np.exp(-1j*self.omega[-1]*t)[None,:,None]
        return values.real

    def initial_bounds(self, radius, depth, ds, dz, deadline):
        active = np.flatnonzero(self.data['carrier_amplitudes_m'])
        C = self.C[active].reshape(len(active),-1)
        k = np.r_[self.k,-self.k]
        ka = abs(k)
        omega = np.tile(self.omega,2)
        sigma = np.sign(k)
        s = np.linspace(-radius,radius,int(np.ceil(2*radius/ds))+1)
        z = np.linspace(depth,0,int(np.ceil(-depth/dz))+1)
        hs,hz = np.max(np.diff(s)),np.max(np.diff(z))
        E0 = np.exp(1j*k[:,None]*s)
        max_v = np.zeros(len(active)); max_dv = max_v.copy()
        for zz in z:
            if time.monotonic() > deadline:
                raise TimeoutError('Declared train preflight budget exceeded')
            E = E0*np.exp(ka[:,None]*zz)
            uh = ((C*(sigma*omega))@E).real
            uz = ((C*(-1j*omega))@E).real
            d1 = ((C*(1j*ka*omega))@E).real
            d2 = ((C*(k*omega))@E).real
            max_v = np.maximum(max_v,np.sqrt(uh*uh+uz*uz).max(axis=1))
            max_dv = np.maximum(max_dv,np.sqrt(2*(d1*d1+d2*d2)).max(axis=1))
        eta = (C@E0).real
        deta = ((C*(1j*k))@E0).real
        r_s = hs/2
        r_sz = np.hypot(hs,hz)/2
        ehess = np.sum(abs(C)*ka**2,axis=1)
        vhess = 2*np.sqrt(2)*np.sum(abs(C)*omega*ka**2,axis=1)
        eta_bound = abs(eta).max(axis=1)+r_s*abs(deta).max(axis=1)+.5*r_s*r_s*ehess
        v_bound = max_v+r_sz*max_dv+.5*r_sz*r_sz*vhess
        return dict(eta_region_bound_m=float(eta_bound.sum()),velocity_region_bound_m_s=float(v_bound.sum()),
            component_eta_bounds_m=eta_bound.tolist(),component_speed_bounds_m_s=v_bound.tolist(),
            projection_interval_m=[-float(radius),float(radius)],depth_interval_m=[float(depth),0.],
            grid_spacing_s_z_m=[float(hs),float(hz)],profile_grid_shape=[len(s),len(z)],
            bound='sum percomponent interval/rectangle Taylor bounds; triangle inequality; shared spectral matrix products')
