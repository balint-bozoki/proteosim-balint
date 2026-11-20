import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Return the list of cleaved peptides from a protein sequence and cleavage pattern.

    Parameters
    ----------
    protein_seq : str
        Protein sequence to be digested, using one-letter code for amino acids.
    
    cleave_pattern : str
        Regular expression specifying the cleavage pattern.

    min_pep_len : int, optional 
        Minimum length (inclusive) of peptides in the output. Default is 5.
    
    max_pep_len : int, optional
        Maximum length (inclusive) of peptides in the output. Default is 30.

    Returns
    -------
    list 
        List of strings containing the peptide sequences whose length is between min_pep_len and max_pep_len.

    """
    

    peptides = re.split(cleave_pattern, protein_seq)
    filtered = [p for p in peptides if min_pep_len <= len(p) <= max_pep_len]
    
    return filtered


def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Return a dictionary mapping protein ID to list of peptides after digestion using digest_protein_sequence.

    Parameters
    ----------
    protein_map : dict
        Dictionary mapping protein IDs to sequence.
    cleave_pattern : str
        Regular expression specifying the cleavage pattern.
    min_pep_len : int, optional 
        Minimum length (inclusive) of peptides in the output. Default is 5.
    max_pep_len : int, optional
        Maximum length (inclusive) of peptides in the output. Default is 30.

    Returns
    -------
    dict 
        Keys contain the protein IDs, values are lists of strings with the peptide sequences whose lengths are between min_pep_len and max_pep_len. 

    """
    digested_protein_map = {}

    for protein, seq in protein_map.items():
        peptides = digest_protein_sequence(seq, cleave_pattern, min_pep_len=min_pep_len, max_pep_len=max_pep_len)
        digested_protein_map[protein] = peptides

    return digested_protein_map


def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the percentage of positions in protein_seq covered by a list of peptides. 

    Parameters
    ----------
    protein_seq : str
        Protein sequence, using one-letter code for amino acids.
    peptides : list
        List of strings containing the peptide sequences to be mapped onto protein_seq. 
   
    Returns
    -------
    float 
        Percentage (0.0-100.0) of residues in protein_seq covered by peptides.
    """
    covered_res =  [False] * len(protein_seq)                                     # start with each position being not covered
    for peptide in peptides: 
        start = protein_seq.find(peptide)                                         # index of the first occurrence
        while start != -1:
            for pos in range(start, start + len(peptide)):
                covered_res[pos] = True
            start = protein_seq.find(peptide, start+1)                            # new search right after current position (overlaps)  
    score = (sum(covered_res) / len(protein_seq)) * 100
    return score


