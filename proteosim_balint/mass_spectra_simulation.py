import numpy as np
import matplotlib.pyplot as plt

amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Return a dictionary mapping peptide sequence to its molecular mass calculated from residues.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence using one letter code for amino acids.
    amino_acid_mass_dict : dict
        Keys are one letter amino acid codes, values are the corresponding molecular masses.

    Returns
    -------
    dict
        Key stores the peptide sequence, value is its molecular mass. 

    """
    result = {peptide_seq:sum(amino_acid_mass_dict[res] for res in peptide_seq)}

    return result

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Return a dictionary mapping peptide sequences to thei molecular masses using a dictionary of amino acid molecular masses. 
    
    Parameters
    ----------
    peptides : list of str
        List of peptide sequences using one letter code for amino acids.
    amino_acid_mass_dict : dict
        Keys are one letter amino acid codes, values are the corresponding molecular masses.

    Returns
    -------
    dict
        Keys are peptide sequences, values are the corresponsing molecular masses calculated from residues.

    """
    result = {}

    for peptide in peptides: 
        result.update(calculate_mol_mass(peptide, amino_acid_mass_dict))
    
    return result

def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Return a dictionary mapping peptide sequences to their m/z values given a charge value.

    Parameters
    ----------
    peptide_mass_map : dict
        Keys are peptide sequences, values are the corresponding molecular masses. 
    charge : int
        Relative charge of the resulting ions.
    proton_mass : float, optional
        Mass of a proton in Da. Default is 1.007.

    Returns
    -------
    dict
        Keys are peptide sequences, values are their corresponding m/z values.

    """
    result = {}
    for peptide, mass in peptide_mass_map.items():
        result[peptide] = (mass + charge * proton_mass) / charge

    return result

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Construct and display an MS spectrum with given m/z values and random intensities in the specified range.

    Parameters
    ----------
    mz_values : dict
        Values contain the m/z values to be used.
    random_count_range : tup
        Tuple specifying the range from which random instensities should be drawn. Default is (0, 3000).
    seed : int, optional
        Specifies the seed for random intensity value generation. Default is 42.

    Returns
    -------
    None

    """
    np.random.seed(seed)
    low, high = random_count_range
    intensities = np.random.randint(low=low, high=high, size=len(mz_values))
    plt.bar(mz_values.values(), intensities, width=1)
    plt.title('MS spectrum')
    plt.xlabel('m/z value')
    plt.ylabel('simulated ion intensity')
    plt.show()

def fragment_peptide(peptide):
    """
    Return a list containing b- and y ion sequences from a peptide sequence.

    Parameters
    ----------
    peptide : str
        Peptide sequence to be fragmented using one letter code for amino acids. 

    Returns
    -------
    list of str
        Entries are b- and y-ion sequences resulting from fragmentation of peptide.

    """
    result = []
    for i in range(1, len(peptide)):          # b-ions
        result.append(peptide[:i])
    for j in range(-1, -len(peptide), -1):    # y-ions
        result.append(peptide[j:])
    
    return result