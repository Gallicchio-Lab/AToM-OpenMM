import os

curr_dir = os.path.dirname(os.path.abspath(__file__))


def _test_uwham_analysis(tmp_path):
    from atom_openmm.uwham import calculate_uwham

    # fmt: off
    lambda1 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50,
               0.40, 0.30, 0.20, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
    lambda2 = [0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50,
               0.50, 0.50, 0.50, 0.50, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]
    # fmt: on

    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_1")
    ddG, ddG_std, dgbind1, dgbind2, samples = calculate_uwham(
        run_dir, "QB_A02_A09", 70, lambda1=lambda1, lambda2=lambda2
    )

    expected_ddG = -1.0582613105156682
    expected_ddG_std = 0.1967262793670161
    assert abs(ddG - expected_ddG) < 1e-8
    assert abs(ddG_std - expected_ddG_std) < 1e-8
    assert samples == 351
    assert abs(dgbind1 - 16.12549393024899) < 1e-8
    assert abs(dgbind2 - 17.183755240764658) < 1e-8

    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_2")
    ddG, ddG_std, dgbind1, dgbind2, samples = calculate_uwham(
        run_dir, "QB_A02_A09", 70, lambda1=lambda1, lambda2=lambda2
    )

    expected_ddG = -1.9492444343451716
    expected_ddG_std = 0.18548486855089807
    assert abs(ddG - expected_ddG) < 1e-8
    assert abs(ddG_std - expected_ddG_std) < 1e-8
    assert samples == 421
    assert abs(dgbind1 - 15.690178371190054) < 1e-8
    assert abs(dgbind2 - 17.639422805535226) < 1e-8

    run_dir = os.path.join(curr_dir, "TYK2_A02_A09", "TYK2_A02_A09_r0_3")
    ddG, ddG_std, dgbind1, dgbind2, samples = calculate_uwham(
        run_dir, "QB_A02_A09", None, lambda1=lambda1, lambda2=lambda2
    )

    expected_ddG = -1.2625416487680994
    expected_ddG_std = 0.2068400471106031
    assert abs(ddG - expected_ddG) < 1e-8
    assert abs(ddG_std - expected_ddG_std) < 1e-8
    assert samples == 480
    assert abs(dgbind1 - 16.159606073746524) < 1e-8
    assert abs(dgbind2 - 17.422147722514623) < 1e-8
