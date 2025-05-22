# Stefan Doerr, 5/2025 UWHAM python implementation based on the R implementation https://cran.r-project.org/web/packages/UWHAM/index.html
import numpy as np
from scipy.optimize import minimize


def _insert(x, d, x0=0):
    """Insert a value x0 at d-th position of x"""
    if d == 0:
        return np.concatenate(([x0], x))
    else:
        return np.concatenate((x[:d], [x0], x[d:]))


def _obj_fcn(ze, logQ, size, base):
    """
    Objective function for UWHAM optimization

    Args:
        ze: vector of log-normalizing constants (or free energies)
        logQ: log of unnormalized density ratios (over baseline)
        size: individual sample sizes for distributions
        base: baseline index
    """
    N = logQ.shape[1]
    rho = size / N

    ze = _insert(ze, base, 0)

    # Calculate normalized Q values
    Qnorm = np.exp(logQ - ze[:, np.newaxis]) * rho[:, np.newaxis]
    Qsum = np.sum(Qnorm, axis=0)

    val = np.sum(np.log(Qsum)) / N + np.sum(ze * rho)

    # Calculate gradient and hessian
    W = np.delete(Qnorm, base, axis=0).T / Qsum[:, None]
    grad = -np.sum(W, axis=0) / N + np.delete(rho, base)

    O = W.T @ W / N
    hess = -O + np.diag(np.sum(W, axis=0)) / N

    return val, grad, hess


def _uwham(logQ, size=None, label=None, base=None, init=None, fisher=True):
    """
    Unbinned Weighted Histogram Analysis Method

    Args:
        logQ: N x M matrix of log unnormalized densities
        size: vector of length M giving sample sizes for distributions
        label: vector of labels indicating which observation is from which distribution
        base: baseline index (0 to M-1) for distribution whose free energy is set to 0
        init: initial values of free energies
        fisher: if True, variance estimation based on Fisher information
    """
    N, M = logQ.shape

    # Compute size if needed
    if size is None:
        if label is None:
            raise ValueError("Either label or size must be provided")
        size = np.bincount(label, minlength=M)[1:]

    # Check which distributions have samples
    sampled = size > 0

    # Validate inputs
    if N != np.sum(size):
        raise ValueError("Inconsistent sum(size) and logQ.shape[0]")
    if M != len(size):
        raise ValueError("Inconsistent length(size) and logQ.shape[1]")

    # Check size and label consistency
    if label is not None:
        _, counts = np.unique(label[label >= 0], return_counts=True)
        if not np.array_equal(counts, size[sampled]):
            raise ValueError("Inconsistent label and size[sampled]")
    else:
        if fisher is False:
            label = np.repeat(np.where(sampled)[0], size[sampled])
            print(
                "Assume observations are ordered by thermodynamic state if fisher=False and label=None"
            )

    # Set base if not provided
    if base is None:
        base = np.where(sampled)[0][0]
    elif not sampled[base]:
        raise ValueError("Observations from baseline are required")

    # Set initial values if not provided
    if init is None:
        init = np.zeros(M)

    # Get number of sampled distributions
    m = np.sum(sampled)
    rho = size[sampled] / N

    # Calculate baseline index for sampled distributions
    temp = np.zeros(M - 1)
    temp = np.insert(temp, base, 1)
    base0 = np.where(temp[sampled])[0][0]

    # Calculate log density ratios
    logQ = logQ - logQ[:, base : base + 1]
    logQ = logQ.T

    # Optimize using scipy
    result = minimize(
        lambda x: _obj_fcn(x, logQ[sampled], size[sampled], base0)[0],
        init[sampled][np.arange(len(init[sampled])) != base0],
        method="trust-exact",
        jac=lambda x: _obj_fcn(x, logQ[sampled], size[sampled], base0)[1],
        hess=lambda x: _obj_fcn(x, logQ[sampled], size[sampled], base0)[2],
    )

    # Process optimization results
    ze0 = _insert(result.x, base0)
    ze = np.zeros(M)
    ze[sampled] = ze0

    # Calculate normalized weights
    Qnorm = np.exp(logQ - ze[:, None])
    Qsum = np.sum(Qnorm[sampled] * rho[:, None], axis=0)

    W = Qnorm.T / Qsum[:, None]
    z = np.mean(W, axis=0)

    # Check values (should all be 1)
    check = z[sampled]

    # Handle unsampled distributions
    if m < M:
        ze[~sampled] = np.log(z[~sampled])
        W[:, ~sampled] = W[:, ~sampled] / z[~sampled]

    # Return early if no variance estimation needed
    if fisher is None:
        return {
            "ze": ze,
            "W": W,
            "check": check,
            "result": result,
            "size": size,
            "base": base,
        }

    # Variance estimation
    O = W.T @ W / N

    D = np.zeros((M, M))
    D[:, sampled] = O[:, sampled] * rho

    H = D - np.eye(M)
    H = np.delete(np.delete(H, base, 0), base, 1)

    if fisher:
        G = O - D[:, sampled] @ O[sampled]
        G = np.delete(np.delete(G, base, 0), base, 1)

        iHG = -O + np.ones((M, 1)) @ O[base : base + 1]
        iHG = np.delete(np.delete(iHG, base, 0), base, 1)

        Ve = iHG @ np.linalg.solve(H.T, np.eye(H.shape[0])) / N

    else:
        C = np.zeros((m, M))
        for j in range(M):
            C[:, j] = np.bincount(label, weights=W[:, j], minlength=M)[sampled]

        G = O - C.T @ np.diag(rho) @ C
        G = np.delete(np.delete(G, base, 0), base, 1)
        Ve = np.linalg.solve(H, G) @ np.linalg.solve(H.T, np.eye(H.shape[0])) / N

    # Variance vector
    ve = _insert(np.diag(Ve), base)

    # Variance-covariance matrix
    Ve = np.insert(Ve, base, 0, axis=1)
    Ve = np.insert(Ve, base, 0, axis=0)

    return {
        "ze": ze,
        "ve": ve,
        "Ve": Ve,
        "W": W,
        "check": check,
        "result": result,
        "label": label,
        "size": size,
        "base": base,
    }


def _bias_fcn(epert, lam1, lam2, alpha, u0, w0):
    """
    Bias ilogistic potential:
    (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    """
    ebias1 = np.zeros_like(epert)
    if alpha > 0:
        ee = 1 + np.exp(-alpha * (epert - u0))
        ebias1 = (lam2 - lam1) * np.log(ee) / alpha
    return ebias1 + lam2 * epert + w0


def _dbias_fcn(epert, lam1, lam2, alpha, u0, w0):
    """
    Derivative of the bias ilogistic potential:
    (lambda2-lambda1)/(1+exp(-alpha (u-u0))) + lambda2
    """
    ee = 1.0 + np.exp(-alpha * (epert - u0))
    return (lam2 - lam1) / ee + lam1


def _npot_fcn(e0, epert, bet, lam1, lam2, alpha, u0, w0):
    """
    Negative reduced energy: -beta*(U0+bias)
    """
    return -bet * (e0 + _bias_fcn(epert, lam1, lam2, alpha, u0, w0))


def _uwham_r(label, logQ, ufactormax, ufactormin=1):
    """
    UWHAM implementation
    Note: This is a simplified version and would need a proper UWHAM implementation
    """
    n = logQ.shape[0]
    m = logQ.shape[1]
    iniz = np.zeros(m)
    uf = ufactormax

    while uf >= ufactormin and uf >= 1:
        mask = slice(0, n, int(uf))
        out = _uwham(label=label[mask], logQ=logQ[mask], init=iniz)
        iniz = out["ze"]
        uf = uf / 2

    return out


def calculate_uwham(
    rundir,
    jobname,
    mintimeid=None,
    maxtimeid=None,
    intermd=None,
    lambda1=None,
    lambda2=None,
    alpha=None,
    u0=None,
    w0=None,
):
    import pandas as pd
    import os

    # Define states
    tempt = np.array([300])
    bet = 1.0 / (0.001986209 * tempt)
    # fmt: off
    # directn = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
    if intermd is None:
        intermd = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    else:
        intermd = np.array(intermd)
    nstates = len(intermd)
    if lambda1 is None:
        lambda1 = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.50, 
                            0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00])
    else:
        lambda1 = np.array(lambda1)
    if lambda2 is None:
        lambda2 = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.50, 
                            0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00])
    else:
        lambda2 = np.array(lambda2)
    if alpha is None:
        alpha = np.full(nstates, 0.10)
    else:
        alpha = np.array(alpha)
    if u0 is None:
        u0 = np.full(nstates, 110.0)
    else:
        u0 = np.array(u0)
    if w0 is None:
        w0 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    else:
        w0 = np.array(w0)
    # fmt: on

    # Calculate states
    # np.where returns a tuple, we want the first element, then take the first occurrence
    leg1istate = np.where(intermd == 1)[0][0]
    leg2istate = np.where(intermd == 1)[0][1]

    # Define column names
    columns = [
        "stateid",
        "temperature",
        "direction",
        "lambda1",
        "lambda2",
        "alpha",
        "u0",
        "w0",
        "potE",
        "pertE",
        "trash",
    ]

    # Create list of data files
    datafiles = [
        os.path.join(rundir, f"r{i}", f"{jobname}.out") for i in range(nstates)
    ]

    # Read and combine all data files
    dfs = []
    for i, file in enumerate(datafiles):
        # Read data file
        df = pd.read_csv(file, sep=r"\s+", header=None, names=columns, index_col=False)
        # Add timeid column
        df["timeid"] = np.arange(1, len(df) + 1)
        dfs.append(df)

    # Combine all dataframes
    data = pd.concat(dfs, ignore_index=True)

    # Add beta column
    data["bet"] = 1.0 / (0.001986209 * data["temperature"])

    # Calculate time masks
    timemask = np.ones(len(data), dtype=bool)
    if mintimeid is not None:
        timemask = timemask & (data["timeid"] >= mintimeid)
    if maxtimeid is not None:
        timemask = timemask & (data["timeid"] <= maxtimeid)

    if not np.any(timemask):
        raise ValueError(
            f"No data found in the given time range. Time values range "
            f"between {data['timeid'].min()} and {data['timeid'].max()} "
            f"and the user requested data between {mintimeid} and {maxtimeid}."
        )

    # Calculate samples
    nsamples = len(data[timemask])
    samplesperreplica = nsamples // nstates

    # Filter data for leg1
    data1 = data[timemask & (data["stateid"] <= leg1istate)]

    mtempt = len(bet)
    leg1stateids = np.arange(leg1istate + 1)
    mlam = len(leg1stateids)
    m = mlam * mtempt
    N = len(data1)

    # Extract U0 values as U-bias
    e0 = data1["potE"].values.copy()
    for i in range(N):
        e0[i] -= _bias_fcn(
            data1["pertE"].iloc[i],
            data1["lambda1"].iloc[i],
            data1["lambda2"].iloc[i],
            data1["alpha"].iloc[i],
            data1["u0"].iloc[i],
            data1["w0"].iloc[i],
        )

    # Calculate negative potential
    neg_pot = np.zeros((N, m))
    sid = 0
    for be in leg1stateids:
        for te in range(mtempt):
            neg_pot[:, sid] = _npot_fcn(
                e0,
                data1["pertE"].values,
                bet[te],
                lambda1[be],
                lambda2[be],
                alpha[be],
                u0[be],
                w0[be],
            )
            sid += 1

    # the alchemical state indexes start with 0, UWHAM's state labels start with 1
    statelabels = data1["stateid"].values.astype(int) + 1

    # Run UWHAM (placeholder - needs implementation)
    out = _uwham_r(label=statelabels, logQ=neg_pot, ufactormax=1, ufactormin=1)

    # Reshape results
    ze = np.array(out["ze"]).reshape(mtempt, mlam)
    ve = np.array(out["ve"]).reshape(mtempt, mlam)

    # Calculate binding free energies
    dgbind1 = (-ze[:, -1] / bet) - (-ze[:, 0] / bet)
    ddgbind1 = np.sqrt(ve[:, -1] + ve[:, 0]) / bet

    # Filter data for leg 2
    data1 = data[timemask & (data["stateid"] >= leg2istate)]

    # Create reversed sequence for leg2 states
    leg2stateids = np.arange(leg2istate, nstates)[::-1]

    mlam = len(leg2stateids)
    m = mlam * mtempt
    N = len(data1)

    # Extract U0 values as U-bias
    e0 = data1["potE"].copy()
    for i in range(N):
        e0.iloc[i] -= _bias_fcn(
            data1["pertE"].iloc[i],
            data1["lambda1"].iloc[i],
            data1["lambda2"].iloc[i],
            data1["alpha"].iloc[i],
            data1["u0"].iloc[i],
            data1["w0"].iloc[i],
        )

    # Calculate negative potential matrix
    neg_pot = np.zeros((N, m))
    sid = 0
    for be in leg2stateids:
        for te in range(mtempt):
            neg_pot[:, sid] = _npot_fcn(
                e0=e0,
                epert=data1["pertE"],
                bet=bet[te],
                lam1=lambda1[be],
                lam2=lambda2[be],
                alpha=alpha[be],
                u0=u0[be],
                w0=w0[be],
            )
            sid += 1

    # Calculate state labels (reversed for leg2)
    statelabels = nstates - data1["stateid"]

    # Run UWHAM
    out = _uwham_r(label=statelabels, logQ=neg_pot, ufactormax=1, ufactormin=1)

    # Reshape results
    ze = out["ze"].reshape(mtempt, mlam)
    ve = out["ve"].reshape(mtempt, mlam)

    # Calculate binding free energies
    dgbind2 = (-ze[:, -1] / bet) - (-ze[:, 0] / bet)
    ddgbind2 = np.sqrt(ve[:, -1] + ve[:, 0]) / bet

    dgb = dgbind1 - dgbind2
    ddgb = np.sqrt(ddgbind2 * ddgbind2 + ddgbind1 * ddgbind1)

    return dgb[0], ddgb[0], dgbind1[0], dgbind2[0], samplesperreplica
