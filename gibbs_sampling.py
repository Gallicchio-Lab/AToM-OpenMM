"""Gibbs sampling routines"""
from numpy import zeros, exp, sum, log, asarray
from numpy.random import random as _random
from random import choice
from itertools import permutations

def _exit(message):
    """Print and flush a message to stdout and then exit."""
    print message
    sys.stdout.flush()
    print 'exiting...'
    sys.exit(1)

def weighted_choice(choices):
    """Return a discrete outcome given a set of outcome/weight pairs."""
    r = _random()*sum(w for c,w in choices)
    for c,w in choices:
        r -= w
        if r < 0:
            return c
    # You should never get here.
    return None

def pairwise_metropolis_sampling(repl_i, sid_i, replicas, states, U):
    """
    Return a replica "j" to exchange with the given replica "i" based on
    the Metropolis criterion:

    P(i<->j) = min{1, exp(-du_ij)};

    du_ij = u_a(x_j) + u_b(x_i) - u_a(x_i) - u_b(x_j),

    where i and j currently occupy states a and b respectively. Repeating this
    MANY times (perhaps N^3 - N^5, N = # of replicas) will construct a 
    Markov chain that will eventually approach the same distribution obtained
    by directly sampling the distribution of all replica/state permutations.
    """
    # Choose another replica other than repl_i.
    #
    nreplicas = len(replicas)
    repl_j = repl_i
    while repl_j == repl_i:
        j = choice(range(nreplicas))
        repl_j = replicas[j]
        sid_j = states[j]
    # Apply the Metropolis acceptance criteria. If the move is accepted, return
    # this replica, otherwise return the same replica (no exchange).
    #
    du = (U[sid_i][repl_j] + U[sid_j][repl_i]
          - U[sid_i][repl_i] - U[sid_j][repl_j])
    if du > 0.:
        if _random() > exp(-du):
            return repl_i
        else:
            return repl_j
    else:
        return repl_j
    
def pairwise_independence_sampling(repl_i, sid_i, replicas, states, U):
    """
    Return a replica "j" to exchange with the given replica "i" based on
    independent sampling from the discrete Metropolis transition matrix, T:

    T_rs = alpha_rs min[1,exp(-du_rs)]  r != s
    T_rr = 1 - sum_(s != r) T_rs        otherwise

    Here r and s are state index permutations, r being the current permutation
    and s the new permutation. alpha_rs = 0 unless r and s differ by a single 
    replica swap and alpha_rs = 1/(n-1) otherwise, n being the number of 
    replicas/states and (n-1) is the number of permutations, s, differing by 
    permutation r by a single swap. du_rs is the change in reduced potential 
    energy of the replica exchange ensemble in going from permutation r to 
    permutation s (that is, due to a replica swap). Based on the above we have:

    du_rs = u_a(j)+u_b(i)-[u_a(i)+u_b(j)]

    where i and j are the replicas being swapped and a and b are, respectively,
    the states they occupy in r and b and a, respectively, those in s.

    The energies u_a(i), i=1,n and a=1,n, are assumed stored in the input 
    "swap matrix," U[a][i].

    In general, the set of replicas across which exchanges are considered is a 
    subset of the n replicas. This list is passed in the 'replicas' list. 
    Replica "i" ('repl_i') is assumed to be in this list.
    """
    # Evaluate all i-j swap probabilities.
    #
    nreplicas = len(replicas)
    ps = zeros(nreplicas) # probability of swap i <-> j
    du = zeros(nreplicas) # Boltzmann exponent, ps ~ exp(-du)
    for j,repl_j,sid_j in zip(range(nreplicas),replicas,states):
        du[j] = (U[sid_i][repl_j] + U[sid_j][repl_i] 
                 - U[sid_i][repl_i] - U[sid_j][repl_j])
    eu = exp(-du)
    
    pii = 1.0
    i = -1
    f = 1./(float(nreplicas) - 1.)
    for j in range(nreplicas):
        repl_j = replicas[j]
        if repl_j == repl_i:
            i = j
        else:
            if eu[j] > 1.0:
                ps[j] = f
            else:                    
                ps[j] = f*eu[j]
            pii -= ps[j]
    try:
        ps[i] = pii
    except IndexError:
        _exit('gibbs_re_j(): unrecoverable error: replica %d not in the '
              'list of waiting replicas?'%i)
    return replicas[weighted_choice(zip(range(nreplicas),ps))]

def state_perm_distribution(replicas, states, swap_matrix):
    """
    Return the distribution of state permutations of a set of replicas (and 
    their states) given the "swap matrix" containing the energies of of all
    such permutations:

    p(s) = (1/Z_s) exp(-sum_i^states u_si)
    
    Z_s = sum_r^permutations exp(-sum_i^states u_ri)

    s/r denote permutations of states, i denotes a specific state
    """
    perm_dist = {}
    Z_s = 0.
    for state_perm in permutations(states):
        perm = str(zip(replicas,state_perm))
        u_s = sum([swap_matrix[j][i] for i,j in zip(replicas,state_perm)])
        exp_us = exp(-u_s)
        perm_dist[perm] = exp_us
        Z_s += exp_us
    for perm,exp_us in perm_dist.iteritems():
        perm_dist[perm] /= Z_s
    return perm_dist

def sample_to_state_perm_distribution(sample, replicas, states):
    """
    Convert a sample count of state permutations to a normalized probability
    distribution of states (including unobserved states).
    """
    Z_s = 0.
    full_sample = {}
    for state_perm in permutations(states):
        perm = str(zip(replicas,state_perm))
        if sample.has_key(perm):
            full_sample[perm] = sample[perm]
        else:
            full_sample[perm] = 0
        Z_s += full_sample[perm]
    for perm,ps in full_sample.iteritems():
        full_sample[perm] /= Z_s
    return full_sample

def state_perm_divergence(empirical, exact):
    """
    Return the Kullback-Liebler divergence of two dicts containing 
    replica/state permutations and the observed probability.
    """
    p = []
    q = []
    for perm,ps_exact in exact.iteritems():
        p.append(empirical[perm])
        q.append(ps_exact)
    return discrete_Kullback_Liebler_divergence(p,q)

def discrete_Kullback_Liebler_divergence(p, q, dx=1.0, eps=1e-9):
    """
    Return the Kullback-Liebler divergence of two discrete, one-dimensional
    probability distributions.

    For numerical stability, change very small elements of the observed
    distribution to some small value, eps. This will yield very large values
    instead of nan or inf due to division by zero.
    """
    p = asarray(p)
    q = asarray(q)
    p[p < eps] = eps
    return (p*log(p/q)*dx).sum()
