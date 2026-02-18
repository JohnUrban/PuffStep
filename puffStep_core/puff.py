import numpy as np
from scipy.stats import norm, poisson, expon, geom, gamma

# Utility functions
def impute_zeros(x, y, bw):
    from scipy.stats import gaussian_kde
    k = gaussian_kde(x, bw_method=bw)
    y_new = np.copy(y)
    y_new[y == 0] = k(x)[y == 0]
    return y_new

def mednorm(x):
    return x / np.median(x)

def mednorm_ksmooth(x, y, bw):
    from scipy.stats import gaussian_kde
    k = gaussian_kde(x, bw_method=bw)
    return mednorm(k(y))

def mednorm_ksmooth_norm(x, y, bw, norm_y):
    return mednorm_ksmooth(x, y, bw) / mednorm_ksmooth(x, norm_y, bw)

def inner_quant_mean_norm(x, inner=(0.4, 0.6)):
    innerq = np.quantile(x, inner)
    return x / np.mean(x[(x >= innerq[0]) & (x <= innerq[1])])

def slope(x1, x2, y1, y2):
    return (y2 - y1) / (x2 - x1)

def prob_data(Forward):
    return np.sum(Forward[:, -1])

# HMM Functions
# Viterbi for normal, poisson, exponential, geometric, gamma, discrete

def viterbi_normal(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + norm.logpdf(emitted_data[0], loc=emissions[0, :], scale=emissions[1, :])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = norm.logpdf(emitted_data[j], loc=emissions[0, i], scale=emissions[1, i]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def viterbi_poisson(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + poisson.logpmf(np.round(emitted_data[0]), mu=emissions[0, :])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = poisson.logpmf(np.round(emitted_data[j]), mu=emissions[0, i]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def viterbi_exponential(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + expon.logpdf(emitted_data[0], scale=1/emissions[0, :])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = expon.logpdf(emitted_data[j], scale=1/emissions[0, i]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def viterbi_geometric(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + geom.logpmf(np.round(emitted_data[0]), p=emissions[0, :])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = geom.logpmf(np.round(emitted_data[j]), p=emissions[0, i]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def viterbi_gamma(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + gamma.logpdf(emitted_data[0], a=emissions[0, :], scale=emissions[1, :])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = gamma.logpdf(emitted_data[j], a=emissions[0, i], scale=emissions[1, i]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def viterbi_discrete(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[0], len(emitted_data)
    pointer = np.zeros((num_emits, num_states), dtype=int)
    Viterbi = np.zeros((num_states, num_emits))
    Viterbi[:, 0] = np.log(initial) + np.log(emissions[:, int(emitted_data[0])])
    pointer[0, :] = 0
    for j in range(1, num_emits):
        selection = Viterbi[:, j-1][:, None] + np.log(transitions)
        for i in range(num_states):
            maxstate = np.argmax(selection[:, i])
            Viterbi[i, j] = np.log(emissions[i, int(emitted_data[j])]) + selection[maxstate, i]
            pointer[j, i] = maxstate
    return Viterbi, pointer

def traceback_viterbi(Viterbi, pointer):
    """Traceback the Viterbi path from the Viterbi matrix and pointer."""
    num_emits = Viterbi.shape[1]
    viterbi_path = np.zeros(num_emits, dtype=int)
    viterbi_path[num_emits-1] = np.argmax(Viterbi[:, num_emits-1])
    
    for j in range(num_emits-1, 0, -1):
        viterbi_path[j-1] = pointer[j, viterbi_path[j]]
    
    # Return 1-based indexing
    return viterbi_path + 1

def posterior(Forward, Backward):
    """Compute posterior decoding from forward and backward matrices."""
    num_emits = Forward.shape[1]
    posterior_path = np.zeros(num_emits, dtype=int)
    for i in range(num_emits):
        fb = Forward[:, i] * Backward[:, i]
        posterior_path[i] = np.argmax(fb) + 1  # +1 for 1-based indexing
    return posterior_path

def forward_normal(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    Forward = np.zeros((num_states, num_emits))
    scalefactors = np.zeros((2, num_emits))
    
    # Initial
    Forward[:, 0] = initial * norm.pdf(emitted_data[0], loc=emissions[0, :], scale=emissions[1, :])
    scalefactors[0, 0] = np.sum(Forward[:, 0])
    scalefactors[1, 0] = np.log(scalefactors[0, 0])
    Forward[:, 0] = Forward[:, 0] / scalefactors[0, 0]
    
    # Iterate
    for k in range(1, num_emits):
        emit = norm.pdf(emitted_data[k], loc=emissions[0, :], scale=emissions[1, :])
        Forward[:, k] = emit * (Forward[:, k-1] @ transitions)
        scalefactors[0, k] = np.sum(Forward[:, k])
        scalefactors[1, k] = np.log(scalefactors[0, k]) + scalefactors[1, k-1]
        Forward[:, k] = Forward[:, k] / scalefactors[0, k]
    
    return Forward, scalefactors

def forward_discrete(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[0], len(emitted_data)
    Forward = np.zeros((num_states, num_emits))
    scalefactors = np.zeros((2, num_emits))
    
    # Initial
    Forward[:, 0] = initial * emissions[:, int(emitted_data[0])]
    scalefactors[0, 0] = np.sum(Forward[:, 0])
    scalefactors[1, 0] = np.log(scalefactors[0, 0])
    Forward[:, 0] = Forward[:, 0] / scalefactors[0, 0]
    
    # Iterate
    for k in range(1, num_emits):
        emit = emissions[:, int(emitted_data[k])]
        Forward[:, k] = emit * (Forward[:, k-1] @ transitions)
        scalefactors[0, k] = np.sum(Forward[:, k])
        scalefactors[1, k] = np.log(scalefactors[0, k]) + scalefactors[1, k-1]
        Forward[:, k] = Forward[:, k] / scalefactors[0, k]
    
    return Forward, scalefactors

def backward_normal(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[1], len(emitted_data)
    Backward = np.zeros((num_states, num_emits))
    scalefactors = np.zeros((2, num_emits))
    
    # Initial
    Backward[:, num_emits-1] = 1
    scalefactors[0, num_emits-1] = np.sum(Backward[:, num_emits-1])
    scalefactors[1, num_emits-1] = np.log(scalefactors[0, num_emits-1])
    Backward[:, num_emits-1] = Backward[:, num_emits-1] / scalefactors[0, num_emits-1]
    
    # Iterate
    for k in range(num_emits-2, -1, -1):
        emit = norm.pdf(emitted_data[k+1], loc=emissions[0, :], scale=emissions[1, :])
        Backward[:, k] = transitions @ (Backward[:, k+1] * emit)
        scalefactors[0, k] = np.sum(Backward[:, k])
        scalefactors[1, k] = np.log(scalefactors[0, k]) + scalefactors[1, k+1]
        Backward[:, k] = Backward[:, k] / scalefactors[0, k]
    
    return Backward, scalefactors

def backward_discrete(emissions, transitions, initial, emitted_data):
    num_states, num_emits = emissions.shape[0], len(emitted_data)
    Backward = np.zeros((num_states, num_emits))
    scalefactors = np.zeros((2, num_emits))
    
    # Initial
    Backward[:, num_emits-1] = 1
    scalefactors[0, num_emits-1] = np.sum(Backward[:, num_emits-1])
    scalefactors[1, num_emits-1] = np.log(scalefactors[0, num_emits-1])
    Backward[:, num_emits-1] = Backward[:, num_emits-1] / scalefactors[0, num_emits-1]
    
    # Iterate
    for k in range(num_emits-2, -1, -1):
        emit = emissions[:, int(emitted_data[k+1])]
        Backward[:, k] = transitions @ (Backward[:, k+1] * emit)
        scalefactors[0, k] = np.sum(Backward[:, k])
        scalefactors[1, k] = np.log(scalefactors[0, k]) + scalefactors[1, k+1]
        Backward[:, k] = Backward[:, k] / scalefactors[0, k]
    
    return Backward, scalefactors
