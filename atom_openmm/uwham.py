# Stefan Doerr, 5/2025 UWHAM python implementation based on the R implementation https://cran.r-project.org/web/packages/UWHAM/index.html
import os
import yaml
import argparse
import glob
import numpy as np
from scipy.optimize import minimize
import pandas as pd
import matplotlib.pyplot as plt

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

def _bias_fcn(epert, par):
    """
    Bias ilogistic potential:
    (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    or
    Multi softplus function    
    depending on par composition
    """
    lam1 = par['lambda1']
    lam2 = par['lambda2']
    alpha = par['alpha']
    u0 = par['u0']
    w0 = par['w0']
    
    if 'lambda3' in par.keys():
        lam3 =  par['lambda3']
        u1 = par['u1']
        f3 = lam3*epert + w0
        if alpha > 0:
            f2 = lam2*epert + w0 + (lam3 - lam2)*u1
            f1 = lam1*epert + w0 + (lam3 - lam2)*u1 + (lam2 - lam1)*u0
            nn = (np.exp(alpha*f1) + np.exp(alpha*f2) + np.exp(alpha*f3))/3.0
            return np.log(nn)/alpha
        else:
            return f3
    else:
        ebias1 = np.zeros_like(epert)
        if alpha > 0:
            ee = 1 + np.exp(-alpha * (epert - u0))
            ebias1 = (lam2 - lam1) * np.log(ee) / alpha
        return ebias1 + lam2 * epert + w0


def _dbias_fcn(epert, par):
    """
    Derivative of the bias ilogistic potential:
    (lambda2-lambda1)/(1+exp(-alpha (u-u0))) + lambda2
    or the Multi softplus function    
    depending on par composition
    """
    lam1 = par['lambda1']
    lam2 = par['lambda2']
    alpha = par['alpha']
    u0 = par['u0']
    w0 = par['w0']
    if 'lambda3' in par.keys():
        lam3 =  par['lambda3']
        u1 = par['u1']
        if alpha > 0:
            f3 = lam3*epert + w0
            f2 = lam2*epert + w0 + (lam3 - lam2)*u1
            f1 = lam1*epert + w0 + (lam3 - lam2)*u1 + (lam2 - lam1)*u0
            nn = (np.exp(alpha*f1) + np.exp(alpha*f2) + np.exp(alpha*f3))/3.0
            return ((lam1*np.exp(alpha*f1) + lam2*np.exp(alpha*f2) + lam3*np.exp(alpha*f3))/3.0)/nn
        else:
            return lam3*np.ones_like(epert)
    else:
        ee = 1.0 + np.exp(-alpha * (epert - u0))
        return (lam2 - lam1) / ee + lam1


def _npot_fcn(e0, epert, bet, par):
    """
    Negative reduced energy: -beta*(U0+bias)
    """
    return -bet * (e0 + _bias_fcn(epert, par))

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

def get_alchemical_schedule(df):
    # create a working copy for deduplication logic
    temp = df.copy()
    
    # deduplicate
    unique_df = temp.drop_duplicates(subset=['stateid'])
    
    # sort by stateid
    unique_df = unique_df.sort_values('stateid')

    # build the alchemical schedule
    if 'lambda3' in unique_df.keys() and 'u1' in unique_df.keys():
        schedule = {
            'stateid': unique_df['stateid'].tolist(),
            'direction': unique_df['direction'].tolist(),
            'lambda1': unique_df['lambda1'].tolist(),
            'lambda2': unique_df['lambda2'].tolist(),
            'lambda3': unique_df['lambda3'].tolist(),
            'alpha':   unique_df['alpha'].tolist(),
            'u0':      unique_df['u0'].tolist(),
            'u1':      unique_df['u1'].tolist(),
            'w0':      unique_df['w0'].tolist(),
            'temperature': unique_df['temperature'].tolist(),
        }        
    else:
        schedule = {
            'stateid': unique_df['stateid'].tolist(),
            'direction': unique_df['direction'].tolist(),
            'lambda1': unique_df['lambda1'].tolist(),
            'lambda2': unique_df['lambda2'].tolist(),
            'alpha':   unique_df['alpha'].tolist(),
            'u0':      unique_df['u0'].tolist(),
            'w0':      unique_df['w0'].tolist(),
            'temperature': unique_df['temperature'].tolist(),
        }

    return schedule

    
def calculate_uwham_leg_from_dataframe(dataf):

    schedule = get_alchemical_schedule(dataf)
    stateid = schedule['stateid']
    vpar = {}
    for lvar in [ "lambda1", "lambda2", "alpha", "u0", "w0"]:
        vpar[lvar] = schedule[lvar]
    do_multi = False
    if 'lambda3' in schedule.keys():
        do_multi = True
        for lvar in [ "lambda3", "u1" ]:
            vpar[lvar] = schedule[lvar]
    else:
        vpar['lambda3'] = schedule['lambda2']
        vpar['u1'] = schedule['u0']
    temperature = schedule['temperature']

    #check that there is only one direction
    unique_direction = list(set(schedule['direction']))
    assert len(unique_direction) == 1, "Each leg data must have only one direction"
    direction = unique_direction[0]
    
    #check that there is only one temperature
    tol = 1.e-9
    temperature_np = np.array(temperature)
    assert np.all(np.abs(temperature_np - temperature_np[0]) < tol), "Multiple temperatures are not currently supported"
    
    mtempt = 1 #only one temperature is supported
    mlam = len(schedule['stateid'])
    m = mlam * mtempt
    N = len(dataf)

    # Add beta column
    kb = 0.001986209
    dataf["bet"] = 1.0 / (kb * dataf["temperature"])
    bet = [ 1.0 / (kb * temperature[0] ) ]
    
    if direction < 0:
        #list states in reverse order for leg2
        stateid = stateid[::-1]
        for lvar in vpar.keys():
            vpar[lvar] = vpar[lvar][::-1]
        temperature = temperature[::-1]

    # Extract U0 values as U-bias
    e0 = dataf["potE"].values.copy()
    for i in range(N):
        par = {}
        for lvar in [ "lambda1", "lambda2", "alpha", "u0", "w0"]:
            par[lvar] = dataf[lvar].iloc[i]
        if do_multi:
            for lvar in [ "lambda3", "u1" ]:
                par[lvar] = dataf[lvar].iloc[i]
        e0[i] -= _bias_fcn( dataf["pertE"].iloc[i], par)

    # Calculate negative potential
    neg_pot = np.zeros((N, m))
    sid = 0
    for be in range(mlam):
        for te in range(mtempt):
            par = {}
            for lvar in [ "lambda1", "lambda2", "alpha", "u0", "w0"]:
                par[lvar] = vpar[lvar][be]
            if do_multi:
                for lvar in [ "lambda3", "u1" ]:
                    par[lvar] = vpar[lvar][be]
            neg_pot[:, sid] = _npot_fcn(
                e0,
                dataf["pertE"].values,
                bet[te],
                par)
            sid += 1

    # assign state labels from 1 to M, with 1 corresponding to lambda=0 or lambda=1 depending on direction
    if direction > 0:
        # the alchemical state indexes start with 0, UWHAM's state labels start with 1
        statelabels = dataf["stateid"].values.astype(int) + 1
    else:
        base = int(np.max(dataf["stateid"].values.astype(int)))
        statelabels = (base - dataf["stateid"].values.astype(int)) + 1
        
    # Run UWHAM (placeholder - needs implementation)
    out = _uwham_r(label=statelabels, logQ=neg_pot, ufactormax=1, ufactormin=1)

    # Reshape results
    ze = np.array(out["ze"]).reshape(mtempt, mlam)
    ve = np.array(out["ve"]).reshape(mtempt, mlam)

    # Calculate binding free energies
    dg = (-ze[:, -1] / bet) - (-ze[:, 0] / bet)
    ddg = np.sqrt(ve[:, -1] + ve[:, 0]) / bet

    return dg[0], ddg[0], out

def get_pertE_leg_distributions(dataf, nbins = 100):
    # A dictionary to hold your results
    densities = {}

    for sid in dataf['stateid'].unique():
        # Filter data for this state
        data = dataf.loc[dataf['stateid'] == sid, 'pertE'].to_numpy()
    
        # Calculate density directly
        # density=True automatically handles the (count / (total_count * bin_width)) calculation
        hist, bin_edges = np.histogram(data, bins=nbins, density=True)
    
        # Calculate bin centers for easier plotting or integration later
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
        densities[sid] = {'centers': bin_centers, 'density': hist}

    return densities

def get_lnp0u_lambdaf_leg(dataf, h = None, n_grid = 100):
    u_grid = np.linspace(dataf['pertE'].min(), dataf['pertE'].max(), n_grid)
    du = u_grid[1] - u_grid[0] 
    
    if h is None:
        h = du

    u_samples = dataf['pertE'].to_numpy()
    weights = dataf['W'].to_numpy()
    #assumes one temperature
    temperature = dataf.iloc[0]['temperature']
    kb = 0.001986209
    bet = 1.0 / (kb * temperature )
    
    # Broadcast subtraction: (len(u_grid), len(u_samples))
    diff = u_grid[:, np.newaxis] - u_samples[np.newaxis, :]
    
    # Gaussian kernel K(x)
    k = (1.0 / np.sqrt(2 * np.pi * h**2)) * np.exp(-0.5 * (diff / h)**2)
    
    # Probability density p(u)
    p = np.sum(weights * k, axis=1) / np.sum(weights)
    
    # Derivative of density p'(u)
    # p'(u) = sum(W_i * K'(x) * (1/h)) where K'(x) = -x * K(x)
    p_prime = np.sum(weights * k * (-diff / h**2), axis=1) / np.sum(weights)
    
    # Log-derivative (lambda function)
    dlnp_du = p_prime / (bet*p)
    
    return u_grid, np.log(p), dlnp_du

def create_quality_assessment_plot(df1, df2):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    prop_cycle = plt.rcParams['axes.prop_cycle']
    default_colors = prop_cycle.by_key()['color']
    nc = len(default_colors)

    plu1 = get_pertE_leg_distributions(df1, nbins = 20)
    plu2 = get_pertE_leg_distributions(df2, nbins = 20)
        
    axs[0, 0].set_title("Leg1 Plambda(u) densities")
    axs[0, 0].set_ylabel("Prob. density (1/(kcal/mol))")
    for i, sid in enumerate(sorted(plu1.keys())):
        axs[0, 0].plot(plu1[sid]['centers'], plu1[sid]['density'], color=default_colors[i % nc], label=f"l{sid}")
    axs[0, 1].set_title("Leg2 Plambda(u) densities")
    for i, sid in enumerate(sorted(plu2.keys(), reverse=True)):
        axs[0, 1].plot(plu2[sid]['centers'], plu2[sid]['density'], color=default_colors[i % nc], label=f"l{sid}")
        
    u_grid_leg1, logp0u_leg1, lambdaf_leg1 = get_lnp0u_lambdaf_leg(df1)
    u_grid_leg2, logp0u_leg2, lambdaf_leg2 = get_lnp0u_lambdaf_leg(df2)

    schedule1 = get_alchemical_schedule(df1)
    schedule2 = get_alchemical_schedule(df2)

    do_multi = False
    if 'lambda3' in schedule1.keys():
        do_multi = True
        
    axs[1, 0].set_title("Leg1 Lambda function")
    axs[1, 0].set_ylim(-0.1, 1.0)
    axs[1, 0].set_ylabel("kBT*dln[p0(u)]/du")
    axs[1, 0].set_xlabel("u (kcal/mol)")
    axs[1, 0].plot(u_grid_leg1, lambdaf_leg1, color='black')
    stateid = schedule1['stateid']

    vpar = {}
    for lvar in [ "lambda1", "lambda2", "alpha", "u0", "w0"]:
        vpar[lvar] = schedule1[lvar]
    if do_multi:
        for lvar in [ "lambda3", "u1" ]:
            vpar[lvar] = schedule1[lvar]
    
    for i in range(len(stateid)):
        par = {}
        for lvar in vpar.keys():
            par[lvar] = vpar[lvar][i]
        axs[1, 0].plot(u_grid_leg1, _dbias_fcn(u_grid_leg1, par), color=default_colors[i % nc])

    axs[1, 1].set_title("Leg2 Lambda function")
    axs[1, 1].set_ylim(-0.1, 1.0)
    axs[1, 1].set_xlabel("u (kcal/mol)")
    axs[1, 1].plot(u_grid_leg2, lambdaf_leg2, color='black')
    stateid = schedule2['stateid'][::-1]

    vpar = {}
    for lvar in [ "lambda1", "lambda2", "alpha", "u0", "w0"]:
        vpar[lvar] = schedule2[lvar][::-1]
    if do_multi:
        for lvar in [ "lambda3", "u1" ]:
            vpar[lvar] = schedule2[lvar][::-1]

    for i in range(len(stateid)):
        par = {}
        for lvar in vpar.keys():
            par[lvar] = vpar[lvar][i]
        axs[1, 1].plot(u_grid_leg2, _dbias_fcn(u_grid_leg2, par), color=default_colors[i % nc])

    return fig
        
def calculate_uwham_from_rundir(
    rundir,
    jobname,
    mintimeid=None,
    maxtimeid=None,
):

    #count number of replicas = number of states looking for folders named r0, r1, etc.
    pattern = os.path.join(rundir, 'r[0-9]*')
    matches = glob.glob(pattern)
    directories = [d for d in matches if os.path.isdir(d)]
    nstates = len(directories)

    #measure the width of the data
    with open(f'{directories[0]}/{jobname}.out', 'r') as f:
        first_line = f.readline()
        ncolumns = len(first_line.split())

    # 11 columns = standard softplus
    # 13 columns = 3-linear softplus
    assert ncolumns == 11 or ncolumns == 13, f"Invalid number of columns in data file: {ncolumns}"
        
    if ncolumns == 11: #softplus
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
            "biasE",
        ]
    else:
        columns = [
            "stateid",
            "temperature",
            "direction",
            "lambda1",
            "lambda2",
            "lambda3",
            "alpha",
            "u0",
            "u1",
            "w0",
            "potE",
            "pertE",
            "biasE",
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

    #extract alchemical schedule
    schedule = get_alchemical_schedule(data)

    # find the stateid's corresponding to the alchemical intermediates
    leg1istate = leg2istate = None 
    ltol = 1.e-6
    lval = 0.5
    inter1 = [ sid for sid in schedule['stateid'] if abs(schedule['lambda1'][sid] - lval) < ltol and abs(schedule['lambda2'][sid] - lval) < ltol and schedule['direction'][sid] == 1 ]
    assert len(inter1) > 0, "Could not find the leg1 alchemical intermediate" 
    assert len(inter1) == 1, "Found multiple leg1 alchemical intermediates"
    leg1istate = inter1[0]
    inter2 = [ sid for sid in schedule['stateid'] if abs(schedule['lambda1'][sid] - lval) < ltol and abs(schedule['lambda2'][sid] - lval) < ltol and schedule['direction'][sid] == -1 ]
    assert len(inter2) > 0, "Could not find the leg2 alchemical intermediate" 
    assert len(inter2) == 1, "Found multiple leg2 alchemical intermediates"
    leg2istate = inter2[0]
    
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
    data1 = data[timemask & (data["stateid"] <= leg1istate)].copy()

    #free energy for leg1
    dg1, ddg1, uwham_out1 = calculate_uwham_leg_from_dataframe(data1)

    # Filter data for leg 2
    data2 = data[timemask & (data["stateid"] >= leg2istate)].copy()

    #free energy for leg2
    dg2, ddg2, uwham_out2 = calculate_uwham_leg_from_dataframe(data2)

    dgb = dg1 - dg2
    ddgb = np.sqrt(ddg1 * ddg1 + ddg2 * ddg2)

    data = {
        'dg_leg1': dg1,
        'dg_stderr_leg1': ddg1,
        'dg_leg2': dg2,
        'dg_stderr_leg2': ddg2,
        'nsamples': samplesperreplica,
        'df_leg1': data1,
        'uwham_out_leg1': uwham_out1,
        'df_leg2': data2,
        'uwham_out_leg2': uwham_out2,
        'schedule': schedule
    }

    return dgb, ddgb, data

def main():
    parser = argparse.ArgumentParser()

    # required input
    parser.add_argument('--jobname', type=str,
                        help='The basename of the job to analyze.', required = True)

    # optional input
    parser.add_argument('--rundir', type=str,
                        help='Processes the data files in the specified directory')
    parser.add_argument('--mintimeid', type=int,
                        help='Process only data obtained after this time id')
    parser.add_argument('--maxtimeid', type=int,
                        help='Process only data obtained before this time id')
    parser.add_argument('--leg1DataCSVoutFile', type=str,
                        help='Saves the samples for leg1 that were processed in a CSV file, includes the WHAM weights')
    parser.add_argument('--leg2DataCSVoutFile', type=str,
                        help='Saves the samples for leg2 that were processed in a CSV file, includes the WHAM weights')
    parser.add_argument('--plotOutFile', type=str,
                        help='A png file where to save the plot for simulation quality assessment')
    args = vars(parser.parse_args())

    jobname = args['jobname']
    
    if args['rundir']:
        rundir = args['rundir']
    else:
        rundir = os.environ.get('PWD')

    discard = None
    if not args['mintimeid']:
        #guess the amount do discard as 1/3 of the data
        outfile0 = os.path.join(rundir, f"r0", f"{jobname}.out")
        with open(outfile0, 'r') as f:
            count = len(f.readlines())
        discard = int(count / 3)
    else:
        discard = args['mintimeid']

    dg, dg_stderr, uwham_data = \
        calculate_uwham_from_rundir(rundir, jobname,
                                    mintimeid=discard, maxtimeid = args['maxtimeid'])
    #save perturbation energy data
    df1 = uwham_data['df_leg1']
    N = len(uwham_data['uwham_out_leg1']['W'][:,0])
    df1['W'] = uwham_data['uwham_out_leg1']['W'][:,0]/float(N)
    if args['leg1DataCSVoutFile']:
        df1.to_csv(args['leg1DataCSVoutFile'], index=False)

    df2 = uwham_data['df_leg2']
    N = len(uwham_data['uwham_out_leg1']['W'][:,0])
    df2['W'] = uwham_data['uwham_out_leg2']['W'][:,0]/float(N)
    if args['leg2DataCSVoutFile']:
        df2.to_csv(args['leg2DataCSVoutFile'], index=False)

    #produces a plot for simulation quality assessment
    if args['plotOutFile']:
        fig = create_quality_assessment_plot(df1, df2)
        fig.savefig(args['plotOutFile'])

    dg1 = uwham_data['dg_leg1']
    dg1_stderr = uwham_data['dg_stderr_leg1']
    dg2 = uwham_data['dg_leg2']
    dg2_stderr = uwham_data['dg_stderr_leg2']
    samplesperreplica = uwham_data['nsamples']
    print(f"{jobname}: DG = {dg:8.3f} +/- {dg_stderr:8.3f}  DG(leg1) = {dg1:8.3f} +/- {dg1_stderr:8.3f}   DG(leg2) = {dg2:8.3f} +/- {dg2_stderr:8.3f} kcal/mol, n_samples: {samplesperreplica}")
    
if __name__ == "__main__":
    main()
