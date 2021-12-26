from stable_fluids import *
import math

"""Stable fluids solver"""
class Simulation:
    """Input:
        dt: time step of the fluid solver.
        n, dissipation, diffusion, viscosity: to initialize a fluid with grid."""
    def __init__(self, dt, n, dissipation, diffusion, viscosity):
        self.dt = dt
        self.fluid = Fluid(n, dissipation, diffusion, viscosity)

    """Reset the simulation, simply reset the fluid to initial state"""
    def reset(self):
        self.fluid = Fluid(self.fluid.n, self.fluid.a, self.fluid.diff, self.fluid.visc)

    def simulate(self):
        self.sim_velocity()
        self.sim_density()

    """Solve for density field"""
    def sim_density(self):
        self.diffusion(self.fluid.d0, self.fluid.d, self.fluid.diff, 0)
        self.advection(self.fluid.d, self.fluid.d0, self.fluid.Vx, self.fluid.Vy, 0)
        self.dissipation()

    """Solve for velocity field"""
    def sim_velocity(self):
        self.diffusion(self.fluid.Vx0, self.fluid.Vx, self.fluid.visc, 1)
        self.diffusion(self.fluid.Vy0, self.fluid.Vy, self.fluid.visc, 2)

        self.projection(self.fluid.Vx0, self.fluid.Vy0, self.fluid.Vx, self.fluid.Vy)

        self.advection(self.fluid.Vx, self.fluid.Vx0, self.fluid.Vx0, self.fluid.Vy0, 1)
        self.advection(self.fluid.Vy, self.fluid.Vy0, self.fluid.Vx0, self.fluid.Vy0, 2)

        self.projection(self.fluid.Vx, self.fluid.Vy, self.fluid.Vx0, self.fluid.Vy0)

    """scalar dissipation, in our case is the density dissipation.
    The dissipation rate is defined as 1+dt*a"""
    def dissipation(self):
        for i in range(self.fluid.s * self.fluid.s):
            self.fluid.d[i] = self.fluid.d[i] / (1 + self.dt * self.fluid.a)

    """Diffusion step of the Navier-Stokes equation.
    The grid will defuse with the diffusion rate dt*k*n*n,
    """
    def diffusion(self, x, x0, k, b):
        diff_rate = self.dt * k * self.fluid.n * self.fluid.n
        self.gauss_seidel(x, x0, diff_rate, 1 + 4 * diff_rate, b)

    """Foster and Metaxas explicit time step with the relaxation schemes.
    They reported Good results with very small iteration of relaxation steps."""
    def gauss_seidel(self, x, x0, a, frac, b):
        for k in range(4):
            for i in range(1, self.fluid.n + 1):
                for j in range(1, self.fluid.n + 1):
                    x[self.fluid.ind(i, j)] = (x0[self.fluid.ind(i, j)] + a * (x[self.fluid.ind(i - 1, j)] +
                                                                               x[self.fluid.ind(i + 1, j)] +
                                                                               x[self.fluid.ind(i, j - 1)] +
                                                                               x[self.fluid.ind(i, j + 1)])) / frac
            self.bnd_cond(x, b)

    """for each grid cell, trace backward in the velocity field, 
    then linear interpolate it and assign back to the grid cell."""
    def advection(self, d, d0, u, v, b):
        for i in range(1, self.fluid.n + 1):
            for j in range(1, self.fluid.n + 1):
                x = i - self.dt * self.fluid.n * u[self.fluid.ind(i, j)]
                y = j - self.dt * self.fluid.n * v[self.fluid.ind(i, j)]

                if x < 0.5:
                    x = 0.5
                if x > self.fluid.n + 0.5:
                    x = self.fluid.n + 0.5
                i0 = math.floor(x)
                i1 = i0 + 1

                if y < 0.5:
                    y = 0.5
                if y > self.fluid.n + 0.5:
                    y = self.fluid.n + 0.5
                j0 = math.floor(y)
                j1 = j0 + 1

                s1 = x - i0
                s0 = 1 - s1
                t1 = y - j0
                t0 = 1 - t1

                d[self.fluid.ind(i, j)] = s0 * (
                        t0 * d0[self.fluid.ind(i0, j0)] + t1 * d0[self.fluid.ind(i0, j1)]) + s1 * (
                                                  t0 * d0[self.fluid.ind(i1, j0)] + t1 * d0[self.fluid.ind(i1, j1)])
        self.bnd_cond(d, b)

    """Projection ensures that velocity is mass conserving.
    Note that this is achieved by using Hodge decomposition. w = u + grad(p)
    and that it can be solved as u = Pw = w - grad(p)
    First we use FDM to calculate the divergence field of velocity,
    then linear solve the Poisson equation and return a incompressible field u where grad(u) = 0.
    last, to obtain u, we have u = w - grad(p) = velocity field - gradient field"""
    def projection(self, u, v, p, div):
        for i in range(1, self.fluid.n + 1):
            for j in range(1, self.fluid.n + 1):
                div[self.fluid.ind(i, j)] = -0.5 * (
                            u[self.fluid.ind(i + 1, j)] - u[self.fluid.ind(i - 1, j)] +
                            v[self.fluid.ind(i, j + 1)] - v[self.fluid.ind(i, j - 1)]) / self.fluid.n
                p[self.fluid.ind(i, j)] = 0
        self.bnd_cond(div, 0)
        self.bnd_cond(p, 0)
        self.gauss_seidel(p, div, 1, 4, 0)

        for i in range(1, self.fluid.n + 1):
            for j in range(1, self.fluid.n + 1):
                u[self.fluid.ind(i, j)] -= 0.5 * (p[self.fluid.ind(i + 1, j)] - p[self.fluid.ind(i - 1, j)]) * self.fluid.n
                v[self.fluid.ind(i, j)] -= 0.5 * (p[self.fluid.ind(i, j + 1)] - p[self.fluid.ind(i, j - 1)]) * self.fluid.n
        self.bnd_cond(u, 1)
        self.bnd_cond(v, 2)

    """Fixed boundary condition for the fluid simulation.
    x is the field we want to fix boundary on, and b is the flag for which boundary condition to fix."""
    def bnd_cond(self, x, b):
        for i in range(1, self.fluid.n + 1):
            if b == 1:
                x[self.fluid.ind(0, i)] = -x[self.fluid.ind(1, i)]
                x[self.fluid.ind(self.fluid.n + 1, i)] = -x[self.fluid.ind(self.fluid.n, i)]
            else:
                x[self.fluid.ind(0, i)] = x[self.fluid.ind(1, i)]
                x[self.fluid.ind(self.fluid.n + 1, i)] = x[self.fluid.ind(self.fluid.n, i)]

            if b == 2:
                x[self.fluid.ind(i, 0)] = -x[self.fluid.ind(i, 1)]
                x[self.fluid.ind(i, self.fluid.n + 1)] = -x[self.fluid.ind(i, self.fluid.n)]
            else:
                x[self.fluid.ind(i, 0)] = x[self.fluid.ind(i, 1)]
                x[self.fluid.ind(i, self.fluid.n + 1)] = x[self.fluid.ind(i, self.fluid.n)]

        x[self.fluid.ind(0, 0)] = 0.5 * x[self.fluid.ind(1, 0)] + x[self.fluid.ind(0, 1)]
        x[self.fluid.ind(0, self.fluid.n + 1)] = 0.5 * x[self.fluid.ind(1, self.fluid.n + 1)] + x[
            self.fluid.ind(0, self.fluid.n)]
        x[self.fluid.ind(self.fluid.n + 1, 0)] = 0.5 * x[self.fluid.ind(self.fluid.n, 0)] + x[
            self.fluid.ind(self.fluid.n + 1, 1)]
        x[self.fluid.ind(self.fluid.n + 1, self.fluid.n + 1)] = 0.5 * x[
            self.fluid.ind(self.fluid.n, self.fluid.n + 1)] + x[self.fluid.ind(self.fluid.n + 1, self.fluid.n)]
