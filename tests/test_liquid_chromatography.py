from proteosim_balint import (predict_lc_retention_times, plot_retention_time, select_retention_time_window)

def test_predict_lc_retention_times():
    peptides = ['PEPTIDE', 'AAA']
    expected = {'PEPTIDE':7.80, 'AAA':6.60}

    actual = predict_lc_retention_times(peptides)
    assert actual == expected

def test_select_retention_time_window():
    peptide_rt_map = {"PEPTIDE":7.80, "AAA":6.60, "AAAPEPTIDE":14.4}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=7, upper_ret_time=10)

    assert selected == ["PEPTIDE"]