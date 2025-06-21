def parse_config(config_file):
    if config_file.endswith(".cntl"):
        from configobj import ConfigObj

        keywords = ConfigObj(config_file, file_error=True)

        # Split float list string keywords
        for key in (
            "DISPLACEMENT",
            "DISPL0_LIG1",
            "DISPL1_LIG1",
            "DISPL0_LIG2",
            "DISPL1_LIG2",
            "LIGOFFSET",
            "TEMPERATURES",
            "LAMBDAS",
            "LAMBDA1",
            "LAMBDA2",
            "ALPHA",
            "U0",
            "W0COEFF",
        ):
            if keywords.get(key, None) is not None:
                keywords[key] = [float(x) for x in keywords[key].split(",")]

        # Split int list string keywords
        for key in ("DIRECTION", "INTERMEDIATE"):
            if keywords.get(key, None) is not None:
                keywords[key] = [int(x) for x in keywords[key].split(",")]

        # Cast list to int
        for key in (
            "LIGAND1_ATOMS",
            "LIGAND2_ATOMS",
            "LIGAND_CM_ATOMS",
            "LIGAND1_CM_ATOMS",
            "LIGAND2_CM_ATOMS",
            "RCPT_CM_ATOMS",
            "POS_RESTRAINED_ATOMS",
            "REST_LIGAND_CMREC_ATOMS",
            "ALIGN_LIGAND1_REF_ATOMS",
            "ALIGN_LIGAND2_REF_ATOMS",
        ):
            if keywords.get(key, None) is not None:
                keywords[key] = [int(x) for x in keywords[key]]

        for key in ("VERBOSE",):
            if keywords.get(key, None) is not None and isinstance(keywords[key], str):
                if keywords[key].lower() == "yes":
                    keywords[key] = True
                elif keywords[key].lower() == "no":
                    keywords[key] = False

        # Try casting the rest to int or float
        for key in keywords:
            if isinstance(keywords[key], str):
                try:
                    keywords[key] = int(keywords[key])
                except ValueError:
                    pass
                try:
                    keywords[key] = float(keywords[key])
                except ValueError:
                    pass

        for key in (
            "MAX_SAMPLES",
            "PRODUCTION_STEPS",
            "PRNT_FREQUENCY",
            "TRJ_FREQUENCY",
            "CHECKPOINT_FREQUENCY",
            "THERMALIZATION_STEPS",
            "ANNEALING_STEPS",
            "EQUILIBRATION_STEPS",
            "STEPS_PER_CYCLE",
            "WALL_TIME",
        ):
            if key in keywords:
                keywords[key] = int(keywords[key])

        return keywords.dict()
    elif config_file.endswith(".yaml") or config_file.endswith(".yml"):
        import yaml

        with open(config_file, "r") as f:
            keywords = yaml.safe_load(f)
    elif config_file.endswith(".json"):
        import json

        with open(config_file, "r") as f:
            keywords = json.load(f)
    else:
        raise ValueError(f"Unsupported config file type: {config_file}")
    return keywords
