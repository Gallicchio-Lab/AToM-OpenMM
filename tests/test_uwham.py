import os

curr_dir = os.path.dirname(os.path.abspath(__file__))


def _test_uwham_analysis(tmp_path):
    from atom_openmm.uwham import calculate_uwham_from_rundir

    dg_tol = 2.e-3 #kcal/mol

    #replicate 1
    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_1")
    ddG, ddG_std, uwham_data = calculate_uwham_from_rundir(
        run_dir, "QB_A02_A09", mintimeid = 70)
    dgbind1 = uwham_data['dg_leg1']
    dgbind2 = uwham_data['dg_leg2']
    samples = uwham_data['nsamples']
    
    expected_n = 351
    expected_ddG = -1.062
    expected_ddG_std = 0.199
    expected_dgbind1 = 15.155
    expected_dgbind2 = 16.2175
    assert samples == expected_n
    assert abs(ddG - expected_ddG) < dg_tol
    assert abs(ddG_std - expected_ddG_std) < dg_tol
    assert abs(dgbind1 - expected_dgbind1) < dg_tol
    assert abs(dgbind2 - expected_dgbind2) < dg_tol

    #replicate 2
    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_2")
    ddG, ddG_std, uwham_data = calculate_uwham_from_rundir(
        run_dir, "QB_A02_A09", mintimeid = 70 )
    dgbind1 = uwham_data['dg_leg1']
    dgbind2 = uwham_data['dg_leg2']
    samples = uwham_data['nsamples']    

    expected_n = 421
    expected_ddG = -1.963
    expected_ddG_std = 0.188
    expected_dgbind1 = 14.709
    expected_dgbind2 = 16.672
    assert samples == expected_n
    assert abs(ddG - expected_ddG) < dg_tol
    assert abs(ddG_std - expected_ddG_std) < dg_tol
    assert abs(dgbind1 - expected_dgbind1) < dg_tol
    assert abs(dgbind2 - expected_dgbind2) < dg_tol

    #replicate 3
    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_3")
    ddG, ddG_std, uwham_data = calculate_uwham_from_rundir(
        run_dir, "QB_A02_A09" )
    dgbind1 = uwham_data['dg_leg1']
    dgbind2 = uwham_data['dg_leg2']
    samples = uwham_data['nsamples']

    expected_n = 480
    expected_ddG = -1.290
    expected_ddG_std = 0.208
    expected_dgbind1 = 15.161
    expected_dgbind2 = 16.451
    assert samples == expected_n
    assert abs(ddG - expected_ddG) < dg_tol
    assert abs(ddG_std - expected_ddG_std) < dg_tol
    assert abs(dgbind1 - expected_dgbind1) < dg_tol
    assert abs(dgbind2 - expected_dgbind2) < dg_tol
