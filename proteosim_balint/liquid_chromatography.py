import matplotlib.pyplot as plt
from pyteomics import achrom


def predict_lc_retention_times(peptides):
    """
    Return a dictionary mapping peptide sequences to relative retention times calculated with achrom.calculate_RT.

    Parameters
    ----------
    peptides: list of str
    List of peptide sequences using one-letter code for amino acids.

    Returns
    -------
    dict
    Keys are peptide sequences, values are the calculated retention times rounded to two decimal places.
    
    """
    results = {}
    for peptide in peptides: 
        results[peptide] = round(achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0), 2)

    return results

def plot_retention_time(retention_times, resolution=30):
    """
    Construct and display a histogram of peptide retention times with the number of bins specified in resolution using matplotlib.

    Parameters
    ----------
    retention_times : list of float
        Peptide retention times in minutes.
    resolutions : int, optional 
        Number of bins in the histogram. Default is 30. 

    Returns
    -------
    None

    """
    plt.hist(retention_times, bins=resolution)
    plt.xlabel("Retention time in minutes")
    plt.ylabel("Number of peptides")
    plt.title("Histogram of predicted retentions times")
    plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Return a list of peptide sequences whose retention times lie between lower_ret_time and upper_ret_time. 

    Parameters
    ----------
    peptide_rt_map : dict
        Keys are peptide sequences, values are retention times in minutes.
    lower_ret_time : float
        Lower threshold (inclusive) of peptide retention times in minutes.
    upper_ret_time : float
        Upper threshold (inclusive) of peptide retention times in minutes.

    Returns
    -------
    list of str
        Entries are peptide sequences whose retention times lie between lower_ret_time and upper_ret_time. 

    """
    filtered = [peptide for peptide, rt in peptide_rt_map.items() if lower_ret_time <= rt <= upper_ret_time]

    return filtered