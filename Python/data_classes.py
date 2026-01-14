import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from constants import *
import my_plot
import dukhin_shilov as ds


class BL_Stat_Data:
    
    def __init__(self, filename=None, Model=None, params=None, path="Data/"):
        
        if params is None and filename is None and Model is None:
            pass
        else:
            self.filename = filename
            self.Model = Model
            self.params = Parameters(params)
            filename = filename + "_bl.txt"
            df = pd.read_csv(path + filename, sep="\s+", skiprows=9, header=None,
                             names=["R", "Z", "U0", "Cp0", "Cn0"])
            df['Distance'] = (df['R'] ** 2 + df['Z'] ** 2) ** 0.5 - self.params.r_c + self.params.r_bl
            df = df.sort_values(by="Distance")
            self.R = df["R"].values
            self.Z = df["Z"].values
            self.U0 = df["U0"].values
            self.Cp0 = df["Cp0"].values
            self.Cn0 = df["Cn0"].values
    
    def plot(self):
        r = np.sqrt(self.R ** 2 + self.Z ** 2)
        fig, ax = plt.subplots(1, 1, figsize=my_plot.set_size("thesis"))
        ax2 = ax.twinx()
        ax.plot(r, self.Cp0, label="$C_{+,0}$")
        ax.plot(r, self.Cn0, label="$C_{-,0}$")
        ax2.plot(r, self.U0, label="$U_0$", color="red")
        ax.set_xlabel("$r [m]$")
        ax.set_ylabel("$C [mol/m^3]$")
        ax.set_yscale("log")
        ax2.set_ylabel("$U_0 [V]$")
        ax2.grid(None)
        # add Brush layer
        ax.axvline(x=self.params.r_c, color="black", linestyle="--")
        ax.axvline(x=self.params.r_c - self.params.r_bl, color="black", linestyle="-")
        from matplotlib.transforms import blended_transform_factory
        x_data_coord = self.params.r_c - self.params.r_bl / 2
        y_axes_coord = 0.5
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(x_data_coord, y_axes_coord, "Brush Layer", transform=trans, rotation=90)
        
        # -- add c_+-, U --
        cp_BL, cn_BL = self.bl_concentrations()
        x = cp_BL / self.params.c0
        U = - kB * T / e * np.log(x)
        x_data_coord = self.params.r_c - self.params.r_bl / 2
        y_data_coord = cp_BL / 1.3
        ax.text(x_data_coord, y_data_coord, "$C_{BL,+}$", va="top", ha="center")
        y_data_coord = cn_BL * 1.3
        ax.text(x_data_coord, y_data_coord, "$C_{BL,-}$", va="bottom", ha="center")
        y_data_coord = U / 1.1
        ax2.text(x_data_coord, y_data_coord, "$U_{BL}$", va="bottom", ha="center")
        
        # add debye length
        c_BL = (cp_BL + cn_BL) / 2 * 0.37
        lambda_D = np.sqrt(eps0 * self.params.epsW * kB * T / (2 * c_BL * NA * e ** 2))
        ax.plot([self.params.r_c, self.params.r_c + lambda_D], [self.params.c0, self.params.c0], color="black",
                lw=0.7)
        ax.text(self.params.r_c + lambda_D / 2, self.params.c0 * 1.2, "$\\lambda_D$", color="black")
        
        fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0.3), bbox_transform=ax.transAxes)
        fig.tight_layout()
        
        return fig, ax, ax2
    
    def to_pdseries(self):
        
        ser = pd.Series()
        ser["filename"] = self.filename
        ser["Model"] = self.Model
        par = self.params.to_pdseries()
        ser = pd.concat([ser, par])
        ser["R"] = self.R
        ser["Z"] = self.Z
        ser["U0"] = self.U0
        ser["Cp0"] = self.Cp0
        ser["Cn0"] = self.Cn0
        
        return ser
    
    def from_pdseries(ser):
        
        data = BL_Stat_Data()
        data.filename = ser["filename"]
        data.Model = ser["Model"]
        data.params = Parameters(ser)
        data.R = ser["R"]
        data.Z = ser["Z"]
        data.U0 = ser["U0"]
        data.Cp0 = ser["Cp0"]
        data.Cn0 = ser["Cn0"]
        
        return data
    
    def __str__(self):
        
        str = f"BL_Stat_Data({self.filename}, {self.params})"
        
        return str
    
    def charge_in_brushlayer(self):
        
        # filter out the brush layer
        R = np.sqrt(self.R ** 2 + self.Z ** 2)
        mask = (R >= self.params.r_c - self.params.r_bl) & (R >= self.params.r_c)
        R = R[mask]
        Cp0 = self.Cp0[mask]
        Cn0 = self.Cn0[mask]
        # calculate the charge in the brush layer
        Q = np.trapz((Cp0 - Cn0) * e * NA, x=R)
        # calculate the charge density
        Sig = Q / (4 * np.pi * (self.params.r_c - self.params.r_bl) ** 2)
        
        return Q / (-self.params.eqSurfCharge)
    
    def bl_el_concentration(self):
        
        R = np.sqrt(self.R ** 2 + self.Z ** 2)
        
        # find closest data point on each side of the interface
        diff = R - self.params.r_c
        diff[diff < 0] = np.inf
        idx_el = np.argmin(diff)
        diff = R - self.params.r_c
        diff[diff > 0] = -np.inf
        idx_bl = np.argmax(diff)
        
        if idx_el == idx_bl:
            idx = np.argmin(np.abs(R - self.params.r_c))
            cp = self.Cp0[idx]
            cn = self.Cn0[idx]
        else:
            def linear_interpolation(x, x1, y1, x2, y2):
                return y1 + ((x - x1) * (y2 - y1)) / (x2 - x1)
            
            cp = linear_interpolation(self.params.r_c, R[idx_el], self.Cp0[idx_el], R[idx_bl], self.Cp0[idx_bl])
            cn = linear_interpolation(self.params.r_c, R[idx_el], self.Cn0[idx_el], R[idx_bl], self.Cn0[idx_bl])
        
        cp_BL, cn_BL = self.bl_concentrations()
        
        c = (cp + cn) / 2
        c_BL = (cp_BL + cn_BL) / 2
        
        exp = False
        if exp:
            cp, cn = (cp-self.params.c0) / (cp_BL-self.params.c0), (cn-self.params.c0) / (cn_BL-self.params.c0)
            c = (c-self.params.c0) / (c_BL-self.params.c0)
        else:
            cp, cn = cp / cp_BL, cn / cn_BL
            c = c / c_BL
        
        return cp, cn, c
    
    def bl_concentrations(self):
        
        cp_BL, cn_BL = self.params.bl_concentrations()
        
        return cp_BL, cn_BL


class Parameters:
    
    def __init__(self, params):
        self.muP = params["muP"]
        self.cP = params["cP"]
        self.epsP = params["epsP"]
        self.mu0 = params["mu0"]
        self.c0 = params["c0"]
        self.epsW = params["epsW"]
        self.r_c = params["r_c"]
        self.L = 10 * self.r_c
        self.r_m = params["r_m"]
        self.epsMem = params["epsMem"]
        self.r_bl = params["r_bl"]
        self.epsBl = params["epsBl"]
        self.muBL = params["muBL"]
        self.eqSurfCharge = params["eqSurfCharge"]
    
    def __str__(self, defaultParams=default_params):
        # check for no-default values
        ser = self.to_pdseries()
        try:
            defaultParams.drop("Model", inplace=True)
        except:
            pass
        ser = ser[ser != defaultParams]
        str = ""
        for key, value in ser.items():
            str += f"{key} = {value:.2e}, "
        
        return str
    
    def to_pdseries(self):
        ser = pd.Series()
        ser["muP"] = self.muP
        ser["cP"] = self.cP
        ser["epsP"] = self.epsP
        ser["mu0"] = self.mu0
        ser["c0"] = self.c0
        ser["epsW"] = self.epsW
        ser["r_c"] = self.r_c
        ser["L"] = self.L
        ser["r_m"] = self.r_m
        ser["epsMem"] = self.epsMem
        ser["r_bl"] = self.r_bl
        ser["epsBl"] = self.epsBl
        ser["muBL"] = self.muBL
        ser["eqSurfCharge"] = self.eqSurfCharge
        
        return ser
    
    def sun_morgan(self, omega, model="Step3"):
        
        # Properties of the surrounding electrolyte
        sig0 = self.mu0 * self.c0 * 2 * F
        epsM_c = eps0 * self.epsW - 1j * sig0 / omega
        lambda_M = np.sqrt(eps0 * self.epsW * kB * T / (2 * self.c0 * NA * e ** 2))
        
        # Properties of the cell plasma
        sigP = self.muP * self.cP * 2 * F
        epsI_c = self.epsP * eps0 - 1j * sigP / omega
        lambda_P = np.sqrt(eps0 * self.epsP * kB * T / (2 * self.cP * NA * e ** 2))
        
        # -- Compute Capacitance of the conceptual membrane --
        
        if model == "Step1":
            cap_P = eps0 * self.epsP / lambda_P  # Capacitance of diffuse layer in cell plasma
            cap_M = eps0 * self.epsW / lambda_M  # Capacitance of diffuse layer in electrolyte
            cap_tot = 1 / (1 / cap_M + 1 / cap_P)  # Total capacitance
            d = lambda_P + lambda_M  # Total thickness of the conceptual membrane
        elif model == "Step2":
            cap_Mem = eps0 * self.epsMem / self.r_m  # Capacitance of the membrane
            cap_P = eps0 * self.epsP / lambda_P  # Capacitance of diffuse layer in cell plasma
            cap_M = eps0 * self.epsW / lambda_M  # Capacitance of diffuse layer in electrolyte
            cap_tot = 1 / (1 / cap_M + 1 / cap_P + 1 / cap_Mem)  # Total capacitance
            d = lambda_P + lambda_M + self.r_m  # Total thickness of the conceptual membrane
        elif model == "Step3":  # modify
            cap_Mem = eps0 * self.epsMem / self.r_m  # Capacitance of the membrane
            cap_P = eps0 * self.epsP / lambda_P  # Capacitance of diffuse layer in cell plasma
            #  Capacitance of diffuse layers at Membrane-BL Interface
            cp_BL, cn_BL = self.bl_concentrations()  # concentrations (deep) in the brush layer
            c_BL = (cp_BL + cn_BL) / 2  # average concentration in the brush layer
            lambda_BL = np.sqrt(eps0 * self.epsP * kB * T / (2 * c_BL * NA * e ** 2))  # Debye length in the brush layer
            cap_BL_Mem = eps0 * self.epsBl / lambda_BL
            #  Capacitance of diffuse layers at BL-EL Interface
            c_BL_EL = c_BL * 0.37  # empirical relation, seems to hold true for C_Bl > 100 C_0
            lambda_BL_EL_inside = np.sqrt(eps0 * self.epsBl * kB * T / (2 * c_BL_EL * NA * e ** 2))  # Debye length at BL-EL Interface
            lambda_BL_EL_outside = np.sqrt(eps0 * self.epsW * kB * T / (2 * c_BL_EL * NA * e ** 2))  # Debye length
            # at BL-EL Interface
            cap_BL_EL_inBL = eps0 * self.epsBl / lambda_BL_EL_inside
            cap_BL_EL_inEL = eps0 * self.epsW / lambda_BL_EL_outside
            cap_tot = 1 / (1 / cap_P + 1 / cap_Mem + 1 / cap_BL_Mem + 1 / cap_BL_EL_inBL + 1 / cap_BL_EL_inEL)
            d = lambda_P + self.r_m + lambda_BL + lambda_BL_EL_inside + lambda_BL_EL_outside
        else:
            raise ValueError("Model not recognized")
        
        # total permitivity of membrane
        eps_Mem = cap_tot * d
        
        # Gamma
        if model == "Step1":
            cellradius = self.r_c
        elif model == "Step2":
            cellradius = self.r_c - self.r_m
        elif model == "Step3":
            cellradius = self.r_c - self.r_bl - self.r_m
        else:
            raise ValueError("Model not recognized")
        
        gamma = (cellradius + d) / cellradius
        
        # permittivity for particle with "Membrane"
        eps_P = eps_Mem * (gamma ** 3 + 2 * (epsI_c - eps_Mem) / (epsI_c + 2 * eps_Mem)) / (
                gamma ** 3 - (epsI_c - eps_Mem) / (epsI_c + 2 * eps_Mem))
        
        # -- Add the surface conductivity of the brush layer --
        if model == "Step3":
            zeta = self.zeta_Potenial(considermuBL=True)  # effective zeta potential
            Sigdp = ds.sigd_vl(self.c0 * NA, self.epsBl, zeta, "+")
            Sigdn = ds.sigd_vl(self.c0 * NA, self.epsBl, zeta, "-")
            kappa = self.mu0 * np.abs(Sigdp - Sigdn)  # surface conductivity of the brush layer
            eps_P = eps_P + 2 * kappa / (1j * omega * cellradius)  # add the surface conductivity to the particle permittivity
        
        # -- Clausius-Mosotti-Factor --
        fcm = (eps_P - epsM_c) / (eps_P + 2 * epsM_c)
        
        # Volume fraction
        phi = 2 * self.r_c ** 3 / (3 * self.L ** 3)
        
        # -- permititvity of mixture --
        epsMix_c = epsM_c * (1 + 2 * phi * fcm) / (1 - phi * fcm)
        
        # -- conductivity of mixture --
        sigMix = epsMix_c * 1j * omega / (self.mu0 * self.c0 * 2 * F + 1j * omega * self.epsW * eps0)
        
        return fcm, sigMix
    
    def zeta_Potenial(self, considermuBL=True):
        
        if considermuBL:
            SigBL = self.charge_in_brushlayer()
            SigEL = 1 - SigBL
            surfCharge = self.eqSurfCharge * (SigBL * self.muBL / self.mu0 + SigEL)
        else:
            surfCharge = self.eqSurfCharge
        
        zeta = 2 * kB * T / e * np.arcsinh(surfCharge / (8 * self.c0 * NA * eps0 * self.epsBl * kB * T) ** 0.5)
        
        return zeta
    
    def maxwell_garnett(self, omega, model="Step3", considermuBL=True):
        
        sigM = 2 * self.c0 * self.mu0 * F + 1j * omega * self.epsW * eps0
        if model == "Step1":
            sigP = 2 * self.cP * self.muP * F + 1j * omega * self.epsP * eps0
        elif model == "Step2":
            sigP = 2 * self.cP * self.muP * F + 1j * omega * self.epsP * eps0
        elif model == "Step3":
            zeta = self.zeta_Potenial(considermuBL)
            Sigdp = ds.sigd_vl(self.c0 * NA, self.epsBl, zeta, "+")
            Sigdn = ds.sigd_vl(self.c0 * NA, self.epsBl, zeta, "-")
            kappa = self.mu0 * np.abs(Sigdp - Sigdn)
            sigP = 2 * self.cP * self.muP * F + 2 * kappa / self.r_c + 1j * omega * self.epsP * eps0
        else:
            raise ValueError("Model not recognized")
        
        nu = 2 * self.r_c ** 3 / (3 * self.L ** 3)
        f = (sigP - sigM) / (sigP + 2 * sigM)
        sigMG = (1 + 2 * nu * f) / (1 - nu * f)
        
        return sigMG
    
    def charge_in_brushlayer(self):
        
        # calculate the ratio of the charges in the brush layer vs the entire diffusive layer using the empirical equation
        old = False
        if old:
            a = 3.13798007e-7
            b = 5.39195671e-3
            x = np.sqrt(self.epsBl / np.abs(self.eqSurfCharge) / self.r_bl)
            Sig_Diff = a * x + b
            Sig_BL = 1 - Sig_Diff
        else:
            a = 3.03934548e-06
            x = np.sqrt(1/np.abs(self.eqSurfCharge) / self.r_bl)
            Sig_Diff = a * x
            Sig_BL = 1 - Sig_Diff
        
        return Sig_BL
    
    def bl_concentrations(self):
        
        rho_bl = self.eqSurfCharge * 3 * (self.r_c - r_bl) ** 2 / (self.r_c ** 3 - (self.r_c - self.r_bl) ** 3)
        c_BL = -rho_bl / (e * NA)  # surplus charge in the brush layer
        
        x = 0.5 * (c_BL / self.c0 + np.sqrt((c_BL / self.c0) ** 2 + 4))
        
        cp_BL = x * self.c0
        cn_BL = 1 / x * self.c0
        
        return cp_BL, cn_BL


class Cond_Spec_Data:
    
    def __init__(self, filename=None, Model=None, params=None, path="Data/"):
        
        if params is None and filename is None and Model is None:
            pass
        else:
            self.filename = filename
            self.Model = Model
            self.params = Parameters(params)
            filename = filename + "_spec.txt"
            df = pd.read_csv(path + filename, sep="\s+", skiprows=5, header=None, names=["omega", "sigma"])
            df['sigma'] = df["sigma"].str.replace('i', 'j').apply(lambda x: complex(x))  # convert to complex Numbers
            self.omega = df["omega"].values
            self.sigma = df["sigma"].values
    
    def __str__(self):
        
        str = f"{self.filename}: {self.params}"
        
        return str
    
    def to_pdseries(self):
        
        ser = pd.Series()
        ser["filename"] = self.filename
        ser["Model"] = self.Model
        par = self.params.to_pdseries()
        ser = pd.concat([ser, par])
        ser["omega"] = self.omega
        ser["sigma"] = self.sigma
        
        return ser
    
    def dukhin_shilov(self, exp=True, omega=None):
        
        sigm = 2 * self.params.c0 * self.params.mu0 * F
        a = self.params.r_c
        mu = self.params.mu0
        C0 = self.params.c0
        epsr = self.params.epsW
        if omega is None:
            w = self.omega
        else:
            w = omega
        
        zeta = self.params.zeta_Potenial(considermuBL=exp)
        
        f_shilov, tau = ds.shilov_vl(sigm, a, w, mu, C0 * NA, epsr, zeta)
        
        nu = 2 * a ** 3 / (3 * self.params.L ** 3)
        
        sigma_shilov = (1 + 2 * nu * f_shilov) / (1 - nu * f_shilov)
        
        return sigma_shilov, tau, zeta
    
    def sun_morgan(self, omega=None):
        
        if omega is None:
            omega = self.omega
        
        fcm, sigMix = self.params.sun_morgan(omega, self.Model)
        
        return fcm, sigMix
    
    def maxwell_garnett(self, considermuBL=True, omega=None):
        
        if omega is None:
            omega = self.omega
        
        sigMG = self.params.maxwell_garnett(omega, self.Model, considermuBL)
        
        return sigMG
    
    def full_ana_solution(self, omega=None):
        
        sigDS, _, _ = self.dukhin_shilov(omega=omega)
        _, sigSM = self.sun_morgan(omega=omega)
        
        sig_Ana = sigDS + sigSM - sigSM[0]
        
        return sig_Ana
    
    def analytical_error(self):
        
        sig_Ana = self.full_ana_solution()
        sig_Num = self.sigma
        
        error_im = np.average(np.abs(np.imag(sig_Ana) - np.imag(sig_Num))) / np.average(np.abs(np.imag(sig_Num)))
        error_re = np.average(np.abs(np.real(sig_Ana) - np.real(sig_Num))) / np.average(np.abs(np.real(sig_Num)))
        
        # correlation coefficient
        corr = np.corrcoef(np.imag(sig_Ana), np.imag(sig_Num))[0, 1]
        
        return error_re, error_im, corr
    
    def rescaled_to_lab(self, nu_lab, sig=None, omega=None, normalize=False):
        
        if (sig is None) and (omega is None):
            sig_0 = self.params.c0 * self.params.mu0 * 2 * F + 1j * self.omega * eps0 * self.params.epsW
            sig = self.sigma * sig_0
        elif (sig is not None) and (omega is not None):
            sig_0 = self.params.c0 * self.params.mu0 * 2 * F + 1j * omega * eps0 * self.params.epsW
            sig = sig * sig_0
        else:
            raise ValueError("Either both sig and omega or none of them must be given")
        
        nu = 2 * self.params.r_c ** 3 / (3 * self.params.L ** 3)
        
        f = 1/nu * (sig - sig_0) / (sig + 2 * sig_0)
        
        sigma_rescaled = sig_0 * (1 + 2 * nu_lab * f) / (1 - nu_lab * f)
        
        if normalize:
            sigma_rescaled = sigma_rescaled / sig_0
        
        return sigma_rescaled
    