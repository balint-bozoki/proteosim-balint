def read_fasta(filepath : str) -> dict:
    """
    Read a FASTA file and return a dictionary mapping protein ID to sequence.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file.
    Returns
    -------
    dict 
        A dictionary mapping protein IDs to sequence.

    """
    protein_map = {}
    current_id = None
    current_sequence = []

    with open(filepath, 'r') as fasta_handle:                                                
        for line in fasta_handle:                                                                
            stripped = line.strip()                                                              
            if not stripped:                                                                     
                continue                                                                         
            if stripped.startswith('>'):                                                         
                if current_id is not None:                                                       
                    protein_map[current_id] = ''.join(current_sequence)                          
                    current_sequence = []
                current_id = stripped.split('|')[1]                                              
            else:
                current_sequence.append(stripped)                                                
    if current_id is not None:                                                                   
        protein_map[current_id] = ''.join(current_sequence)
    return protein_map