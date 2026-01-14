import pandas as pd
import numpy as np

# Fundamental constants
F = 96485.3  # Faraday's constant (C/mol)
kB = 1.3806E-23  # Boltzmann's constant (J/K)
eps0 = 8.854e-12  # vacuum permitivity [A·s/(V·m)]
e = 1.602e-19  # elemtal charge [A*s]
NA = F/e  # Avogadro constant [1/mol]

# Default Model parameters
T = 293  # Absolute temperature [K]
mu = 5e-8  # Ion mobility in bulk electrolyte (m^2/(Vs))
epsW = 80  # Relative permittivity of electrolyte (-)
epsP = 80  # Relative permittivity of solid phase (-)
c0 = 1  # Bulk conc. (mol/m^3)
a = 1e-6  # Particle radius (m)
L = 10 * a  # Radius / half-length of modelling domain (m)
E0 = 1  # External electrical field (V/m)
muP = mu
cP = c0
r_m = 1e-8  # Membrane thickness (m)
epsMem = 5  # Relative permittivity of membrane (-)
r_bl = 1e-8  # Brush layer thickness (m)
epsBl = 80  # Relative permittivity of brush layer (-)
muBL = mu/10  # Ion mobility in brush layer (m^2/(Vs))
eqSurfCharge = -0.5  # equivalent surface charge density (C/m^2)

# Fixed model parameters and derived quantities
D = mu * kB * T  # Ion diffusivity in bulk electrolyte (m^2/s)
ka = 2 * c0 * mu * F  # Electrolyte conductivity (S/m)
nu = 2 * a ** 3 / (3 * L ** 3)  # Volume percent of particles (-)
sigma0 = 2*c0*mu*F  # Conductivity of the electrolyte (S/m)

default_params = pd.Series(dtype=float)
default_params["muP"] = muP
default_params["cP"] = cP
default_params["epsP"] = epsP
default_params["mu0"] = mu
default_params["c0"] = c0
default_params["epsW"] = epsW
default_params["r_c"] = a
default_params["L"] = L
default_params["Model"] = "Step1"
default_params["r_m"] = np.NAN  # Base model does not have a membrane with finite thickness
default_params["epsMem"] = np.NAN  # Base model does not have a membrane with finite thickness
default_params["r_bl"] = np.NAN  # Only used for the brush layer model
default_params["epsBl"] = np.NAN  # Only used for the brush layer model
default_params["muBL"] = np.NAN  # Only used for the brush layer model
default_params["eqSurfCharge"] = np.NAN  # Only used for the brush layer model
