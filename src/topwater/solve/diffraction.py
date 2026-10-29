"""Coherent packet diffraction, shared by motion forcing and water recovery."""
import numpy as np


class PacketDiffraction:
    def __init__(self, hydro, train):
        self.hydro, self.train = hydro, train
        self.results = {}
        self.forces = {}
        theta = np.mod(np.asarray(train.theta)[:, None]+[0., np.pi], 2*np.pi).ravel()
        self.directions = np.column_stack((np.cos(theta), np.sin(theta)))

    def incident(self, index, points, *, velocity=False):
        """Deep-water Airy waves in Capytaine's exp(-i*omega*t) convention."""
        omega = float(self.train.omega[index])
        k = omega**2/self.hydro.gravity
        coefficients = np.asarray(self.train.C[:, :, index]).ravel()
        active = coefficients != 0
        directions, coefficients = self.directions[active], coefficients[active]
        phase = np.exp(k*points[:, 2, None]+1j*k*(points[:, :2]@directions.T))
        if velocity:
            vectors = np.column_stack((directions, np.full(len(directions), -1j)))
            return self.hydro.gravity*k/omega*(phase@(coefficients[:, None]*vectors))
        return self.hydro.density*self.hydro.gravity*(phase@coefficients)

    def solve(self, index):
        if index not in self.results:
            from capytaine.bem.problems_and_results import LinearPotentialFlowProblem
            hydro = self.hydro
            mesh = hydro.body.mesh_including_lid
            velocity = self.incident(index, mesh.faces_centers, velocity=True)
            rhs = -np.einsum('ij,ij->i', velocity, mesh.faces_normals)
            if np.any(rhs):
                problem = LinearPotentialFlowProblem(body=hydro.body,
                    omega=float(self.train.omega[index]), boundary_condition=rhs,
                    rho=hydro.density, g=hydro.gravity, water_depth=hydro.depth)
                self.results[index] = hydro.solver.solve(problem, keep_details=True)
            else:
                self.results[index] = None
        return self.results[index]

    def force(self, index):
        if index not in self.forces:
            from topwater.solve.types import MODES
            result = self.solve(index)
            pressure = self.incident(index, self.hydro.body.mesh.faces_centers)
            fk = self.hydro.body.integrate_pressure(pressure)
            self.forces[index] = np.array([fk[mode]+(0. if result is None else result.forces[mode])
                                          for mode in MODES])
        return self.forces[index]


def packet_diffraction(hydro, train):
    response = getattr(hydro, 'packet_response', None)
    if response is None or response.train is not train:
        response = hydro.packet_response = PacketDiffraction(hydro, train)
    return response
