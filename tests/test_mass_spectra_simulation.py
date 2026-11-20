from proteosim_balint import (calculate_mol_mass, calculate_mol_mass_collection, calculate_mz_collection, fragment_peptide)

def test_calculate_mol_mass():
    aa_mass_dict = {'A':10.5, 'B':20.0, 'C':2.5}
    test_peptide = 'ABC'
    expected_mass = {'ABC':33.0}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

def test_calculate_mol_mass_collection():
    aa_mass_dict = {'A':10.5, 'B':20.0, 'C':2.5}
    peptides = ['AA', 'ABC', 'BBB', 'BCA']
    expected = {'AA':21.0, 'ABC':33.0, 'BBB':60.0, 'BCA':33.0}
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected

def test_calculate_mz_collection():
    peptide_mass_map = {'A':10.0, 'B':20.0, 'C':30.5}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {'A':6.007, 'B':11.007, 'C':16.257}

    assert actual == expected

def test_fragment_peptide():
    peptide = 'PEPT'
    expected = ['P', 'PE', 'PEP', 'T', 'PT', 'EPT']
    actual = fragment_peptide(peptide)

    assert set(actual) == set(expected)
    assert peptide not in expected