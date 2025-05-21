from setuptools import setup

setup(
    scripts=[
        "atom_openmm/abfe_structprep.py",
        "atom_openmm/abfe_production.py",
        "atom_openmm/rbfe_structprep.py",
        "atom_openmm/rbfe_production.py",
        "atom_openmm/make_atm_system_from_pdb.py",
        "atom_openmm/make_atm_system_from_Amber.py",
        "atom_openmm/make_atm_system_from_rcpt_lig.py",
    ]
)
