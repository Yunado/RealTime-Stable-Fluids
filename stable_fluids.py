import numpy as np

"""Fluid class for the stable fluids simulation"""
class Fluid:
    """Input:
        n: size of the grid.
        a: dissipation rate of the fluid.
        diffusion: diffusion constant of the fluid.
        viscosity: kinematic viscosity of the fluid.

        Note that the grid is constructed with 2 extra layers
        for fixed boundary. Thus, the grid has size (n+2)x(n+2)

        Construct a fluid grid where:
        s: Actual size of the grid (n+2)x(n+2)
        Vx: X component of the velocity field.
        Vy: Y component of the velocity field.
        Vx0: Previous step of the X velocity field.
        Vy0: Previous step of the Y velocity field.
        d: the scalr quantity, in our case the density field.
        d0: Previous step of the density field."""
    def __init__(self, n, dissipation, diffusion, viscosity):
        self.n = n
        self.s = n + 2
        self.a = dissipation
        self.diff = diffusion
        self.visc = viscosity
        self.Vx = np.zeros(self.s * self.s)
        self.Vy = np.zeros(self.s * self.s)
        self.Vx0 = np.zeros(self.s * self.s)
        self.Vy0 = np.zeros(self.s * self.s)
        self.d = np.zeros(self.s * self.s)
        self.d0 = np.zeros(self.s * self.s)

    """Get the 1d array index from 2d grid location"""
    def ind(self, x, y):
        return x + y * self.s

    """Breakdown add_force to two component, add density and add velocity"""
    def add_density(self, x, y, amount):
        i = self.ind(x, y)
        self.d[i] += amount

    """Breakdown add_force to two component, add density and add velocity"""
    def add_velocity(self, x, y, amountX, amountY):
        i = self.ind(x, y)
        self.Vx[i] += amountX
        self.Vx[i] += amountY
