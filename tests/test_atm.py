import shutil
import os


curr_dir = os.path.dirname(os.path.abspath(__file__))


def _test_abfe_structprep(tmp_path):
    from atom_openmm.abfe_structprep import (
        massage_keywords,
        do_mintherm,
        do_lambda_annealing,
        do_equil,
    )
    from atom_openmm.utils.config import parse_config
    import logging

    logger = logging.getLogger("abfe_structprep")
    run_dir = os.path.join(tmp_path, "1STP")
    shutil.copytree(os.path.join(curr_dir, "1STP"), run_dir)
    os.chdir(run_dir)

    keywords = parse_config("1STP_asyncre.yaml")
    restrain_solutes = True
    old_keywords = keywords.copy()
    massage_keywords(keywords, restrain_solutes)

    do_mintherm(keywords, logger)
    do_lambda_annealing(keywords, logger)

    # reestablish the restrained atoms
    if restrain_solutes:
        keywords["POS_RESTRAINED_ATOMS"] = old_keywords.get("POS_RESTRAINED_ATOMS")

    do_equil(keywords, logger)


def _test_abfe_production(tmp_path):
    from atom_openmm.openmm_async_re import openmm_job_ABFE

    shutil.copytree(
        os.path.join(curr_dir, "1STP_equil"),
        os.path.join(tmp_path, "1STP_equil"),
    )
    run_dir = os.path.join(tmp_path, "1STP_equil")
    os.chdir(run_dir)

    rx = openmm_job_ABFE("1STP_asyncre.yaml", options=None)
    rx.setupJob()
    rx.scheduleJobs()


def _test_rbfe_structprep(tmp_path):
    from atom_openmm.rbfe_structprep import (
        massage_keywords,
        do_mintherm,
        do_lambda_annealing,
        do_equil,
    )
    from atom_openmm.utils.config import parse_config
    import logging
    from openmm import unit
    from openmm import app
    from openmm.app.amberprmtopfile import AmberPrmtopFile
    from openmm import XmlSerializer
    import numpy as np

    logger = logging.getLogger("rbfe_structprep")
    run_dir = os.path.join(tmp_path, "QB_A08_A07")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07"), run_dir)
    os.chdir(run_dir)

    box_vectors = np.array([[80.1356, 0, 0], [0, 86.0443, 0], [0, 0, 81.4926]])
    a = unit.Quantity(box_vectors[0] * unit.angstrom)
    b = unit.Quantity(box_vectors[1] * unit.angstrom)
    c = unit.Quantity(box_vectors[2] * unit.angstrom)
    box_vectors = (a, b, c)
    omm_structure = AmberPrmtopFile("QB_A08_A07.prmtop", periodicBoxVectors=box_vectors)
    system = omm_structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=9 * unit.angstrom,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,
        rigidWater=True,
        switchDistance=7 * unit.angstrom,
    )
    with open("QB_A08_A07_sys.xml", "w") as f:
        f.write(XmlSerializer.serialize(system))
    keywords = parse_config("QB_A08_A07_asyncre.yaml")

    restrain_solutes = True

    old_keywords = keywords.copy()
    massage_keywords(keywords, restrain_solutes)

    do_mintherm(keywords, logger)
    do_lambda_annealing(keywords, logger)

    # reestablish the restrained atoms
    if restrain_solutes:
        keywords["POS_RESTRAINED_ATOMS"] = old_keywords.get("POS_RESTRAINED_ATOMS")

    do_equil(keywords, logger)


def _test_rbfe_production(tmp_path):
    from atom_openmm.openmm_async_re import openmm_job_RBFE

    shutil.copytree(
        os.path.join(curr_dir, "QB_A08_A07_equil"),
        os.path.join(tmp_path, "QB_A08_A07_equil"),
    )
    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil")
    os.chdir(run_dir)

    rx = openmm_job_RBFE("QB_A08_A07_asyncre.yaml", options=None)
    rx.setupJob()
    rx.scheduleJobs()


def _test_input_parser():
    from atom_openmm.utils.config import parse_config

    config_yaml = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.yaml"))
    config_json = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.json"))
    config_cntl = parse_config(os.path.join(curr_dir, "QB_A08_A07_asyncre.cntl"))

    assert config_yaml == config_json
    assert config_yaml == config_cntl

    for key in config_yaml.keys():
        assert config_yaml[key] == config_json[key]
        assert config_yaml[key] == config_cntl[key]
        assert type(config_yaml[key]) is type(
            config_json[key]
        ), f"{key} {type(config_yaml[key])} {type(config_json[key])}"
        assert type(config_yaml[key]) is type(
            config_cntl[key]
        ), f"{key} {type(config_yaml[key])} {type(config_cntl[key])}"


def _test_sync_production_ats(tmp_path):
    from atom_openmm.openmm_async_re import openmm_job_RBFE

    shutil.copytree(
        os.path.join(curr_dir, "QB_A08_A07_equil_sync_ats"),
        os.path.join(tmp_path, "QB_A08_A07_equil_sync_ats"),
    )

    run_dir = os.path.join(tmp_path, "QB_A08_A07_equil_sync_ats")
    os.chdir(run_dir)

    rx = openmm_job_RBFE("QB_A08_A07_asyncre.yaml", options=None)
    rx.setupJob()
    rx.scheduleJobs()
    for i in range(4):
        assert os.path.exists(
            os.path.join(
                tmp_path, "QB_A08_A07_equil_sync_ats", f"r{i}", "QB_A08_A07.xtc"
            )
        )


def _test_sync_production_incremental(tmp_path):
    from atom_openmm.openmm_async_re import openmm_job_RBFE
    from openmm.app.internal.xtc_utils import read_xtc
    import yaml

    def rbfe_production(configfile):
        rx = openmm_job_RBFE(configfile, options=None)
        rx.setupJob()
        rx.scheduleJobs()

    exec_dir = os.path.join(tmp_path, "QB_A08_A07_equil_sync")
    shutil.copytree(os.path.join(curr_dir, "QB_A08_A07_equil_sync"), exec_dir)

    os.chdir(exec_dir)

    configfile = os.path.join(exec_dir, "QB_A08_A07_asyncre.yaml")
    with open(configfile, "r") as f:
        config = yaml.safe_load(f)
    config["MAX_SAMPLES"] = "+1"
    with open(configfile, "w") as f:
        yaml.dump(config, f)

    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(exec_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 1)
        assert box_read.shape == (3, 3, 1)
        assert len(time) == 1
        assert len(step) == 1

    startsampl_file = os.path.join(exec_dir, "starting_sample")
    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert starting_sample == 0
    prog_file = os.path.join(exec_dir, "progress")
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.0

    # Test restarting
    config["MAX_SAMPLES"] = "+2"
    with open(configfile, "w") as f:
        yaml.dump(config, f)
    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(exec_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 2)
        assert box_read.shape == (3, 3, 2)
        assert len(time) == 2
        assert len(step) == 2

    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert starting_sample == 0  # The starting_sample file existed so was not overwritten
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.5

    # Add two more samples
    os.remove(startsampl_file)
    os.remove(prog_file)
    config["MAX_SAMPLES"] = "+2"
    with open(configfile, "w") as f:
        yaml.dump(config, f)
    rbfe_production(configfile)
    for i in range(4):
        xtcf = os.path.join(exec_dir, f"r{i}", "QB_A08_A07.xtc")
        coords_read, box_read, time, step = read_xtc(xtcf.encode("utf-8"))
        assert coords_read.shape == (46380, 3, 4)
        assert box_read.shape == (3, 3, 4)
        assert len(time) == 4
        assert len(step) == 4

    with open(startsampl_file, "r") as f:
        starting_sample = int(f.read().strip())
        assert starting_sample == 2  # We removed the starting_sample file so it was written from the checkpoint
    with open(prog_file, "r") as f:
        progress = float(f.read().strip())
        assert progress == 0.5
