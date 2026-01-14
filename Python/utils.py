import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit

import cole_cole
from data_classes import BL_Stat_Data, Cond_Spec_Data

from cole_cole import *
from constants import *


def read_saved_data(dict_file, data_dir = "Data/", silent=True):
    """
    Reads the data from the files specified in the dictionary file and returns the data as a list of Data objects.
    :param dict_file: name of the dictionary file (including path)
    :return: Lists of Data objects
    """
    
    # define default parameters
    params = default_params.copy()
    params = params.drop("Model")
    params["cP"] = 10
    params["r_m"] = 1e-8  # Base model does not have a membrane with finite thickness
    params["epsMem"] = 5  # Base model does not have a membrane with finite thickness
    params["r_bl"] = 1e-8  # Only used for the brush layer model
    params["epsBl"] = epsW  # Only used for the brush layer model
    params["muBL"] = mu  # Only used for the brush layer model
    params["eqSurfCharge"] = -0.5  # Only used for the brush layer model
    Model = "Step3"
    
    # read the dictionary file
    df_dict = pd.read_csv(dict_file, sep='\s+', header=0)
    My_Dict = []
    for index, row in df_dict.iterrows():
        Var_params = row["Varied_Parameters"].split(",")
        Var_Values = row["Parameter_Values"].split(",")
        my_dict = dict(zip(Var_params, Var_Values))
        My_Dict.append(my_dict)
    
    Names = df_dict["Filename"]
    BL_Data = []
    Spec_Data = []
    
    for i in range(len(Names)):
        par_temp = params.copy()
        filename = Names[i]
        for key, value in My_Dict[i].items():
            par_temp[key] = value
            par_temp = par_temp.astype(float)
        # check if file exists
        try:
            bl_data = BL_Stat_Data(filename, Model, par_temp, data_dir)
            BL_Data.append(bl_data)
        except:
            if not silent:
                print(f"File {filename}_bl.txt not found")
        try:
            spec_data = Cond_Spec_Data(filename, Model, par_temp, data_dir)
            Spec_Data.append(spec_data)
        except:
            if not silent:
                print(f"File {filename}_spec.txt not found")
    
    return BL_Data, Spec_Data


def get_debye_lengths(c0, cP, T=T, epsW=epsW, epsP=epsP):
    """
    Calculates the Debye lengths of electrolyte and cell plasma
    
    Parameters:
    c0 (float): The bulk concentration of electrolyte.
    cP (float): The bulk concentration of cell plasma.
    
    Returns:
    tuple: The Debye lengths of electrolyte and cell plasma.
    """
    
    lambda_e = np.sqrt(eps0 * epsW * kB * T / (2 * c0 * NA * e ** 2))
    lambda_p = np.sqrt(eps0 * epsP * kB * T / (2 * cP * NA * e ** 2))
    
    return lambda_e, lambda_p


def read_table(name, cols=["w", "var", "s"]):
    """
    Reads a table from a CSV file.

    Parameters:
    name (str): The name of the CSV file.
    cols (list): The column names.

    Returns:
    DataFrame: The data from the CSV file as a pandas DataFrame.
    """
    
    X = pd.read_csv(name, sep="\s+", skiprows=5, header=None, names=cols)
    X['s'] = X['s'].str.replace('i', 'j').apply(lambda x: complex(x))  # convert to complex Numbers
    X["source_file"] = name
    return X


def plot_results(df, labelvar="c", name=None, save=False):
    fig, ax = plt.subplots(1, 1)
    for var in df["var"].unique():
        label = labelvar + " = " + "{:.2e}".format(var)
        ax.scatter(df.loc[df["var"] == var]["w"], np.imag(df.loc[df["var"] == var]["s"]), s=10, label=label)
    ax.set_xscale("log")
    ax.grid("both", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("$\omega$ [Hz]")
    ax.set_ylabel("Im($\sigma$)")
    fig.legend()
    fig.tight_layout()
    if save:        fig.savefig(f"{name}.pdf", bbox_inches='tight')
    
    return fig


def get_peaks(x, type="both"):
    peaks = argrelextrema(x, np.greater_equal, order=2)
    valleys = argrelextrema(x, np.less_equal, order=2)
    if type == "both":
        ind = np.concatenate((peaks, valleys), axis=None)
    elif type == "pos":
        ind = peaks
    elif type == "neg":
        ind = valleys
    ind = np.setdiff1d(ind, [0, len(x) - 1])
    
    return ind


def plot_peaks(df, minx, maxx, type="both"):
    df = df.loc[(df["w"] < maxx) & (df["w"] > minx)]
    Ind = []
    Var = []
    fig, ax = plt.subplots(2, 1)
    for var in df["var"].unique():
        s_im = np.imag(df.loc[df["var"] == var]["s"])
        w = np.array(df.loc[df["var"] == var]["w"])
        ind = get_peaks(s_im, type)
        ax[0].scatter(np.full_like(w[ind], var), w[ind])
        ax[1].scatter(np.full_like(s_im[ind], var), s_im[ind])
        Ind.append(ind)
        Var.append(var)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("$\omega$")
    ax[0].grid("both")
    ax[1].set_ylabel("Im($\sigma$)")
    ax[1].grid("both")
    return Ind, Var, fig


def add_params_to_df(df, params, varying_param=None):
    if varying_param is not None:
        df = df.rename(columns={"var": varying_param})
    
    if not ("Model" in df.columns):
        df["Model"] = params["Model"]
    if not ('muP' in df.columns):
        df["muP"] = params["muP"]
    if not ('cP' in df.columns):
        df["cP"] = params["cP"]
    if not ('epsP' in df.columns):
        df["epsP"] = params["epsP"]
    if not ('mu0' in df.columns):
        df["mu0"] = params["mu0"]
    if not ('c0' in df.columns):
        df["c0"] = params["c0"]
    if not ('epsW' in df.columns):
        df["epsW"] = params["epsW"]
    if not ("r_c" in df.columns):
        df["r_c"] = params["r_c"]
    if not ("L" in df.columns):
        df["L"] = params["L"]
    if not ("r_m" in df.columns):
        df["r_m"] = params["r_m"]
    if not ("epsMem" in df.columns):
        df["epsMem"] = params["epsMem"]
    if not ("r_bl" in df.columns):
        df["r_bl"] = params["r_bl"]
    if not ("epsBl" in df.columns):
        df["epsBl"] = params["epsBl"]
    if not ("muBL" in df.columns):
        df["muBL"] = params["muBL"]
    if not ("eqSurfCharge" in df.columns):
        df["eqSurfCharge"] = params["eqSurfCharge"]
    
    return df


def merge_dataframes(db, df):
    len = db.shape[0]
    
    db = pd.concat([db, df])
    
    param_cols = list(db.columns)
    param_cols.remove("sigma")
    param_cols.remove("omega")
    
    # db = db.sort_values(by=['muP', 'cP', "epsP", "mu0", "c0", "epsW", "r_c", "r_m", "source_file"])
    db = db.sort_values(by=param_cols)
    db.reset_index(drop=True, inplace=True)
    # db.drop_duplicates(['muP', 'cP', "epsP", "mu0", "c0", "epsW", "r_c", "source_file"], inplace = True)
    db.drop_duplicates(param_cols, inplace=True)
    param_cols.remove("source_file")
    # duplicates = db[db.duplicated(['muP', 'cP', "epsP", "mu0", "c0", "epsW", "r_c"], keep = 'first')]
    duplicates = db[db.duplicated(param_cols, keep='first')]
    if duplicates.empty:
        print("No duplicates")
    else:
        print(f"!Double Entries!: {duplicates.shape[0]}")
    print(f"New entries added: {db.shape[0] - len}")
    return db


def collapse_rows(df):
    param_cols = list(df.columns)
    param_cols.remove("s")
    param_cols.remove("w")
    params = df[param_cols].copy()
    params.drop_duplicates(inplace=True)
    
    data = pd.DataFrame(columns=["omega", "sigma"]).astype(object)
    
    for index, row in params.iterrows():
        mask = np.all((df[param_cols] == row) | (df[param_cols].isna()), axis=1)
        data.loc[index, "sigma"] = np.array(df["s"][mask])
        data.loc[index, "omega"] = np.array(df["w"][mask])
    
    params["sigma"] = data["sigma"]
    params["omega"] = data["omega"]
    
    return params


def plot_fit_result(freq, rho, popt, type="cond"):
    fig, axes = plt.subplots(2, 1)
    
    if type == "res":
        axes[0].plot(freq, np.imag(rho), label="data")
        axes[1].plot(freq, np.real(rho))
        axes[0].plot(freq, np.imag(cole_cole.colecole(freq, popt[0], popt[1], popt[2], 1, 1)), label="fit",
                     linestyle="--")
        axes[1].plot(freq, np.real(cole_cole.colecole(freq, popt[0], popt[1], popt[2], 1, 1)), linestyle="--")
        axes[0].set_ylabel("Im($\\rho$)")
        axes[1].set_ylabel("Re($\\rho$)")
    elif type == "cond":
        axes[0].plot(freq, np.imag(1 / rho), label="data")
        axes[1].plot(freq, np.real(1 / rho))
        axes[0].plot(freq, np.imag(1 / cole_cole.colecole(freq, popt[0], popt[1], popt[2], 1, 1)), label="fit",
                     linestyle="--")
        axes[1].plot(freq, np.real(1 / cole_cole.colecole(freq, popt[0], popt[1], popt[2], 1, 1)), linestyle="--")
        axes[0].set_ylabel("Im($\sigma$)")
        axes[1].set_ylabel("Re($\sigma$)")
    
    axes[0].set_xscale("log")
    axes[0].set_xlabel("$\omega$")
    axes[0].grid("both")
    
    axes[1].set_xscale("log")
    axes[1].set_xlabel("$\omega$")
    axes[1].grid("both")
    fig.legend()
    return fig


def filter_df(df, muP=None, cP=None, epsP=None, c0=None, mu0=None, epsW=None):
    """
    Filters a DataFrame based on given conditions. The conditions are checked for approximate equality
    using the numpy.isclose function. The function can handle both list and non-list inputs for each parameter.

    Parameters:
    df (DataFrame): The DataFrame to filter.
    muP, cP, epsP, c0, mu0: The conditions to filter the DataFrame on. Each can be a single value or a list of values.

    Returns:
    DataFrame: The filtered DataFrame.
    """
    if muP is not None:
        if not isinstance(muP, list):
            muP = [muP]
        df = df[np.any([np.isclose(df["muP"], val, rtol=1e-10, atol=0) for val in muP], axis=0)]
    
    if cP is not None:
        if not isinstance(cP, list):
            cP = [cP]
        df = df[np.any([np.isclose(df["cP"], val, rtol=1e-10, atol=0) for val in cP], axis=0)]
    
    if epsP is not None:
        if not isinstance(epsP, list):
            epsP = [epsP]
        df = df[np.any([np.isclose(df["epsP"], val, rtol=1e-10, atol=0) for val in epsP], axis=0)]
    
    if c0 is not None:
        if not isinstance(c0, list):
            c0 = [c0]
        df = df[np.any([np.isclose(df["c0"], val, rtol=1e-10, atol=0) for val in c0], axis=0)]
    
    if mu0 is not None:
        if not isinstance(mu0, list):
            mu0 = [mu0]
        df = df[np.any([np.isclose(df["mu0"], val, rtol=1e-10, atol=0) for val in mu0], axis=0)]
    
    if epsW is not None:
        if not isinstance(epsW, list):
            epsW = [epsW]
        df = df[np.any([np.isclose(df["epsW"], val, rtol=1e-10, atol=0) for val in epsW], axis=0)]
    
    df.reset_index(inplace=True, drop=True)
    
    return df


def get_subcaption_pos(ax):
    """
    Gets the position for a subcaption in a plot.

    Parameters:
    ax (Axes): The axes object to get the subcaption position for.

    Returns:
    tuple: The x and y coordinates for the subcaption.
    """
    bbox = ax.get_position()
    x_middle = bbox.x0 + bbox.width / 2
    y_lower = bbox.y0 - 0.08
    return x_middle, y_lower


def select_simulations(df, params, values, filter_default=True, default_params=default_params, remove_duplicates=True):
    """
    Selects the simulations from a DataFrame that match the given parameter values. Not listed parameters are filtered to the default value.

    Parameters:
    df (DataFrame): The DataFrame to select the simulations from.
    params (list): The parameters to select the simulations on.
    values (list): Values for the parameters to select the simulations on. Can be a single value, a list of values, or a list of lists of values.
    filter_default (bool): Whether to filter the parameters not listed in the selection to the default value.
    default_params (Series): The default parameter values to filter the parameters not listed in the selection to.
    remove_duplicates (bool): Whether to remove duplicate rows by parameter values.

    Returns:
    DataFrame: The selected simulations.
    """
    
    # Function to check if a series is close to any of the given values
    def isclose_to_any(series, values, rtol=1e-10, atol=0):
        """
        Checks if a series is close to any of the given values.
        """
        return np.any([np.isclose(series, val, rtol=rtol, atol=atol) for val in values], axis=0)
    
    # Filter the DataFrame to the given parameter values
    mask = np.all(
        [isclose_to_any(df[param], val) if isinstance(val, list) else np.isclose(df[param], val, rtol=1e-10, atol=0) for
         param, val in zip(params, values)], axis=0)
    df = df[mask]
    df.reset_index(inplace=True, drop=True)
    
    if filter_default:
        # Filter the DataFrame to the default parameter values for the parameters not listed in the selection
        # Get the default parameter names
        def_param_names = default_params.index.tolist()
        for param in params:
            def_param_names.remove(param)
        # filter by string default values
        string_params = [param for param in def_param_names if isinstance(default_params[param], str)]
        mask = np.all([df[param] == df[param] for param in string_params], axis=0)
        df = df[mask]
        df.reset_index(inplace=True, drop=True)
        # remove string parameters from the list
        def_param_names = [param for param in def_param_names if param not in string_params]
        # remove nan values from the list
        def_param_names = [param for param in def_param_names if not np.isnan(default_params[param])]
        # Filter the DataFrame by remaining default parameter values
        for param in def_param_names:
            mask = np.isclose(df[param], default_params[param], rtol=1e-10, atol=0)
            df = df[mask]
            df.reset_index(inplace=True, drop=True)
    
    if remove_duplicates:
        # Remove duplicate rows by parameter values
        df.drop_duplicates(inplace=True, subset=default_params.index.tolist())
    
    return df


def format_number(num):
    num = round(num, 1)
    str_num = str(num)
    formatted_num = str_num.rstrip('0').rstrip('.') if '.' in str_num else str_num
    return formatted_num
