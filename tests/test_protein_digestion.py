from proteosim_balint import (digest_protein_collection, digest_protein_sequence, compute_sequence_coverage, enzyme_cleavage_patterns)

def test_digest_protein_sequence():
    dummy_seq_1 = ""
    dummy_seq_2 = "AAAAAAAAAAAAAAAAAAAAAA"
    dummy_seq_3 = "KPAAAARPAAAKAAAPAAARAA"



    assert digest_protein_sequence(dummy_seq_1, r'(?<=[KR])(?!P)') == []
    assert digest_protein_sequence(dummy_seq_2, r"(?<=K)") == [dummy_seq_2]
    assert digest_protein_sequence(dummy_seq_3, r"(?<=K)") == ["PAAAARPAAAK", "AAAPAAARAA"]
    assert digest_protein_sequence(dummy_seq_3, r"(?<=K)", min_pep_len=1) == ["K", "PAAAARPAAAK", "AAAPAAARAA"]
    assert digest_protein_sequence(dummy_seq_3, r"(?<=K)", max_pep_len=10) == ["AAAPAAARAA"]
    assert digest_protein_sequence(dummy_seq_3, r"(?<=K)", min_pep_len=1, max_pep_len=10) == ["K", "AAAPAAARAA"]
    assert digest_protein_sequence(dummy_seq_3, r'(?<=[KR])(?!P)', min_pep_len=1) == ["KPAAAARPAAAK", "AAAPAAAR", "AA"]

def test_digest_protein_collection():
    dummy_proteins = {
        'PROTEIN_1': 'AAAAAAAAAAAAAAAAAAAAAA', 
        'PROTEIN_2': 'KPAAAARPAAAKAAAPAAARAA', 
        'PROTEIN_3': '', 
        'P11802': 'MATSRYEPVAEIGVGAYGTVYKARDPHSGHFVALKSVRVPNGGGGGGGLPISTVREVALLRRLEAFEHPNVVRLMDVCATSRTDREIKVTLVFEHVDQDLRTYLDKAPPPGLPAETIKDLMRQFLRGLDFLHANCIVHRDLKPENILVTSGGTVKLADFGLARIYSYQMALTPVVVTLWYRAPEVLLQSTYATPVDMWSVGCIFAEMFRRKPLFCGNSEADQLGKIFDLIGLPPEDDWPRDVSLPRGAFPPRGPRPVQSVVPEMEESGAQLLLEMLTFNPHKRISAFRALQHSYLHKDEGNPE'
    }

    assert digest_protein_collection(dummy_proteins, r'(?<=K)') == {'PROTEIN_1': ['AAAAAAAAAAAAAAAAAAAAAA'],
 'PROTEIN_2': ['PAAAARPAAAK', 'AAAPAAARAA'],
 'PROTEIN_3': [],
 'P11802': ['MATSRYEPVAEIGVGAYGTVYK',
  'ARDPHSGHFVALK',
  'VTLVFEHVDQDLRTYLDK',
  'APPPGLPAETIK',
  'DLMRQFLRGLDFLHANCIVHRDLK',
  'PENILVTSGGTVK',
  'PLFCGNSEADQLGK',
  'RISAFRALQHSYLHK',
  'DEGNPE']}
    assert digest_protein_collection(dummy_proteins, r'(?<=[KR])(?!P)', min_pep_len=10, max_pep_len=15) == {'PROTEIN_1': [],
 'PROTEIN_2': ['KPAAAARPAAAK'],
 'PROTEIN_3': [],
 'P11802': ['DPHSGHFVALK',
  'LEAFEHPNVVR',
  'VTLVFEHVDQDLR',
  'APPPGLPAETIK',
  'GLDFLHANCIVHR',
  'KPLFCGNSEADQLGK',
  'IFDLIGLPPEDDWPR']}

def test_compute_sequence_coverage():
    dummy_prot_seq = 'KPRPAAAKWAHAPARTRE'
    dummy_peps = ['WAHAPAR', 'RP']

    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    assert coverage == 50.0 