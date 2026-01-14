from config import *
from utils import *
import utils as ut

# --- read data for Step1+2 ---
df = pd.read_pickle(r"Step_1+2_Data.pkl")  # returns a DataFrame of all conducted simulations of Step 1+2
# containing data and parameters

# --- read data for Step3 (Full model with brush layer) ---
# Additions/Changes to the default parameters
default_params["Model"] = "Step3"
default_params["cP"] = 10
default_params["r_m"] = 1e-8
default_params["epsMem"] = 5
default_params["r_bl"] = 1e-8
default_params["epsBl"] = epsW
default_params["muBL"] = mu
default_params["eqSurfCharge"] = -0.5
# Dictionary with the parameters that are varied
dict_file = "Data_Brushlayer_Dictionary.txt"
# Directory with the data files
data_dir = "Data_Brushlayer/"
BL_Data, Spec_Data = utils.read_saved_data(dict_file, data_dir)  # returns two lists of objects defined in data_classes.py
# BL_Data contains data of a cut through the brush layer (only for a few simulations)
# Spec_Data contains the conductivity spectra of the simulations

# --- calculate analytical solution for Step3 ---
omega = np.logspace(-1, 9, 1000)  # frequency range
Sig_Ana = []  # list to store the analytical solutions
for spec_data in Spec_Data:
    sig_ana = spec_data.full_ana_solution(omega)
    Sig_Ana.append(sig_ana)

