from collections import defaultdict
import numpy as np
from scipy import stats as sps
#from pk2txt import bdgmsg, newmsg
from puffStep_core.utils import *
import puffStep_core.puff as puff

## PREVIOUSLY THIS FUNCTION USED VARIABLES DEFINED INSIDE R.
## NEWER VERSION BELOW USES ONLY VARIABLES DEFINED IN PYTHON FIRST.
##def hmm7(late, path, emodel):
##    states = {}
##    if path == 'viterbi':
##        for chrom in late.chromosomes:
####            sys.stderr.write( chrom + "\n" )
##            if len(late.count[chrom]) > 1:
##                v = puffR.viterbi_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##                states[chrom] = list(v[0])
##            else:
##                states[chrom] = [0] ## if only 1 bin, assign a non-state
##    elif path == 'posterior':
##        for chrom in late.chromosomes:
##            f = puffR.forward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##            b = puffR.backward_puff(emissions = puffR.emissions7, transitions = puffR.transitions7, initial = puffR.initial7, states = intvec([1,2,3,4,5,6,7]), emitted_data = fltvec(late.count[chrom]), emodel = emodel, logprobs=False)
##            p = posterior(f[0], b[0], [1,2,3,4,5,6,7])
##            states[chrom] = list(p[0])            
##    return states


## NEWER VERSION OF HMM7 BELOW USES ONLY VARIABLES DEFINED IN PYTHON -- SEE PROBABILITY MATRICES BELOW.
## IT ALSO ALLOWS SOME TOGGLING OF HMM PARAMETERS

def viterbipath(late, emodel, eprobs, tprobs, iprobs, states):
    statepath = {}
    for chrom in late.chromosomes:
        if len(late.count[chrom]) > 1:
            emitted_data = late.count[chrom]
            if emodel == "normal":
                V, pointer = puff.viterbi_normal(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "discrete":
                V, pointer = puff.viterbi_discrete(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "poisson":
                V, pointer = puff.viterbi_poisson(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "exponential":
                V, pointer = puff.viterbi_exponential(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "geometric":
                V, pointer = puff.viterbi_geometric(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "gamma":
                V, pointer = puff.viterbi_gamma(eprobs, tprobs, iprobs, emitted_data)
            viterbi_path = puff.traceback_viterbi(V, pointer)
            statepath[chrom] = list(viterbi_path)
        else:
            statepath[chrom] = [0]
    return statepath

def posteriorpath(late, emodel, eprobs, tprobs, iprobs, states):
    statepath = {}
    for chrom in late.chromosomes:
        if len(late.count[chrom]) > 1:
            emitted_data = late.count[chrom]
            if emodel == "normal":
                f, _ = puff.forward_normal(eprobs, tprobs, iprobs, emitted_data)
                b, _ = puff.backward_normal(eprobs, tprobs, iprobs, emitted_data)
            elif emodel == "discrete":
                f, _ = puff.forward_discrete(eprobs, tprobs, iprobs, emitted_data)
                b, _ = puff.backward_discrete(eprobs, tprobs, iprobs, emitted_data)
            p = puff.posterior(f, b)
            statepath[chrom] = list(p)
        else:
            statepath[chrom] = [0]
    return statepath

def findstatepath(late, path, emodel, eprobs, tprobs, iprobs, states):
    if path == 'viterbi':
        statepath = viterbipath(late, emodel, eprobs, tprobs, iprobs, states)
    elif path == 'posterior':
        statepath = posteriorpath(late, emodel, eprobs, tprobs, iprobs, states)
    return statepath


def learn_iprobs(late, statepath, nstates, learnpseudo=1e-323):
    inits = [learnpseudo]*nstates
    for chrom in late.chromosomes:
        for state in statepath[chrom]:
            if state > 0:
                inits[int(state)-1] += 1

    np_iprobs = get_initial_probs(inits=inits)
    return np_iprobs

def learn_tprobs(late, statepath, nstates, learnpseudo=1e-323):
    np_tprobs = np.zeros((nstates,nstates)) + learnpseudo
    for chrom in late.chromosomes:
        for x in range(1, len(statepath[chrom])):
            state_i = int(statepath[chrom][x-1])-1
            state_j = int(statepath[chrom][x])-1
            if state_i >= 0 and state_j >= 0:
                np_tprobs[state_i, state_j] += 1
    np_tprobs = (np_tprobs.transpose()/np_tprobs.sum(1)).transpose()
    return np_tprobs

def learn_eprobs(late, statepath, nstates, old_np_eprobs, emodel="normal", learnpseudo=1e-323, constrain=False):
    if not constrain:
        edict = {state:[] for state in range(nstates)}
        for chrom in late.chromosomes:
            for x in range(len(statepath[chrom])):
                state = int(statepath[chrom][x])-1
                if state >= 0:
                    edict[state].append( late.count[chrom][x] )
    if emodel == "discrete":
        nsym = old_np_eprobs.shape[1]
        np_eprobs = old_np_eprobs
        if not constrain:
            for state in range(nstates):
                if edict[state]:
                    edict[state] = np.array(edict[state])
                    symcounts = dict(zip(*np.unique(edict[state], return_counts=True)))
                    for sym in range(nsym):
                        if sym in symcounts.keys():
                            np_eprobs[state, sym] = symcounts[sym]
                        else:
                            np_eprobs[state, sym] = learnpseudo
                else:
                    pass
            np_eprobs = (np_eprobs.transpose()/np_eprobs.sum(1)).transpose()
    else:
        np_eprobs = old_np_eprobs
        if not constrain:
            for state in range(nstates):
                if edict[state]:
                    edict[state] = np.array(edict[state])
                    np_eprobs[0, state] = edict[state].mean()
                    np_eprobs[1, state] = edict[state].std(ddof=1)
                else:
                    pass

    return np_eprobs


def updateparameters(late, statepath, nstates, old_np_eprobs, learnpseudo=1e-323, emodel="normal", emitpseudo=1e-7, constrainEmit=False):
    np_iprobs = learn_iprobs(late, statepath, nstates)
    np_tprobs = learn_tprobs(late, statepath, nstates, learnpseudo=learnpseudo)
    np_eprobs = learn_eprobs(late, statepath, nstates, old_np_eprobs=old_np_eprobs, emodel=emodel, learnpseudo=emitpseudo, constrain=constrainEmit)
    return np_eprobs, np_tprobs, np_iprobs


def log10_prob_state_path_eNormal(late, statepath, np_eprobs, np_tprobs, np_iprobs):
    logprob = 0.0
    for chrom in late.chromosomes:
        ## GET all possible Eprobs in one fell swoop -- believe it or not, this is MUCH faster than getting one at a time in loop
        epd = {es:np.log10(
                sps.norm(np_eprobs[0, es],
                         np_eprobs[1, es]).pdf(
                             late.count[chrom])
                ) for es in range(np_tprobs.shape[0])}
        ## Initial prob and initial emission prob:
        init_state = int(statepath[chrom][0])-1
        if init_state >= 0: # zero reserved for errors (e.g. contigs too small)
            ## add init
            logprob += np.log10(np_iprobs[init_state])
            ## add emit
            ## OLD METH: logprob += np.log10( sps.norm(np_eprobs[0, init_state], np_eprobs[1, init_state]).pdf( late.count[chrom][0]) )
            logprob += epd[init_state][0]
        ## ITER FOR TRANS AND EMIT PROBS
        for x in range(1, len(statepath[chrom])):
            state_i = int(statepath[chrom][x-1])-1
            state_x = int(statepath[chrom][x])-1
            if state_i >= 0 and state_x >= 0:
                ## Add trans
                logprob += np.log10(np_tprobs[state_i, state_x])
                ## Add emit
                ## OLD METH: logprob += np.log10( sps.norm(np_eprobs[0, state_x], np_eprobs[1, state_x]).pdf( late.count[chrom][x]) )
                logprob += epd[state_x][x]
        ## FOR DEVEL: newmsg(chrom + '\t' + str(logprob))
    return logprob

def log10_prob_state_path_eDiscrete(late, statepath, np_eprobs, np_tprobs, np_iprobs):
    logprob = 0.0
    for chrom in late.chromosomes:
        ## Initial prob and initial emission prob:
        init_state = int(statepath[chrom][0])-1
        if init_state >= 0: # zero reserved for errors (e.g. contigs too small)
            ## add init
            logprob += np.log10(np_iprobs[init_state])
            ## add emit
            ## OLD METH: logprob += np.log10( sps.norm(np_eprobs[0, init_state], np_eprobs[1, init_state]).pdf( late.count[chrom][0]) )
            logprob += np.log10( np_eprobs[init_state, 0] )
        ## ITER FOR TRANS AND EMIT PROBS
        for x in range(1, len(statepath[chrom])):
            state_i = int(statepath[chrom][x-1])-1
            state_x = int(statepath[chrom][x])-1
            if state_i >= 0 and state_x >= 0:
                ## Add trans
                logprob += np.log10(np_tprobs[state_i, state_x])
                ## Add emit
                ## "Discrete" assumes the symbols emitted are sequential integers from 1:N where N = number states.
                ## Therefore, for proper pythonese indexing 1 needs to be subtracted from the observed emitted symbol
                emit = int(late.count[chrom][x] - 1) ## 
                logprob += np.log10( np_eprobs[state_x, emit] )
        ## FOR DEVEL: newmsg(chrom + '\t' + str(logprob))
    return logprob

def log10_prob_state_path(late, statepath, np_eprobs, np_tprobs, np_iprobs, emodel):
    if emodel == "normal":
        return log10_prob_state_path_eNormal(late, statepath, np_eprobs, np_tprobs, np_iprobs)
    elif emodel == "discrete":
        return log10_prob_state_path_eDiscrete(late, statepath, np_eprobs, np_tprobs, np_iprobs)


    
def hmmR(late, path, emodel, eprobs=None, tprobs=None, iprobs=None, iters=1):
    if eprobs is None:
        eprobs = emissions7()
    if tprobs is None:
        tprobs = transitions7()
    if iprobs is None:
        iprobs = initial7()
    states = list(range(1, len(iprobs)+1))
    statepath = findstatepath(late, path, emodel, eprobs, tprobs, iprobs, states)           
    return statepath

def generate_hmmR(late, emodel, eprobs=None, tprobs=None, iprobs=None):
    raise NotImplementedError("generate_hmmR is not implemented in pure Python version. R/rpy2 dependency removed.")
    # states = list(range(1, len(iprobs)+1))
    # statepath = {}
    # emitted_data = {}
    # for chrom in late.chromosomes:
    #     statepathlen = len(late.count[chrom])
    #     if statepathlen > 1:
    #         ans = puffR.generate(emissions = eprobs, transitions = tprobs, initial = iprobs, states = states, statepathlen = statepathlen, emodel = emodel)
    #         statepath[chrom] = list(ans[0])
    #         emitted_data[chrom] = list(ans[1])
    #     else:
    #         pass
    # return statepath, emitted_data


## FXNS TO HELP GENERATE PROB MATRICES FOR R
def get_emission_probs(mu=[1,2,4,8,16,32,64],sig=None, discrete=False):
    if discrete:
        np_eprobs = np.array(mu, dtype=float)
    else:
        if sig is None:
            sig = [e**0.5 for e in mu]
        np_eprobs = np.array([mu, sig], dtype=float)
    return np_eprobs


##def get_transition_probs(nstates=7, changestate=0.001, samestate=0.999):
##    t = np.zeros([nstates,nstates], dtype=float)
##    for i in range(nstates):
##      t[i,] = changestate #0.000001
##      t[i,i] = samestate #0.999999
##      t[i,] = t[i,]/sum(t[i,])
##    rowsvec = []
##    for i in range(nstates):
##        rowsvec += list(t[i,])
##    rowsvecr = fltvec(rowsvec)
##    return matrixr(rowsvecr, nrow=nstates, byrow=True)

def get_simple_transition_probs(nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001):
    stay_special = 1 - (leave_special * (nstates-1))
    stay_non = 1 - leave_non_to_special - (leave_non_to_othernon * (nstates-2))
    # nonspecial = list(range(nstates)) ## ERASE LINE???
    # nonspecial.pop(special_state_idx) ## ERASE LINE???
    np_tprobs = np.zeros([nstates,nstates], dtype=float)

    # Matrix is first designed with probabilities reflecting non-special state moves: non-to-othernon, non-to-special, and non-to-selfnon. 
    for i in range(nstates):
        # matrix row i is initialized with leave_non_to_othernon for all other states
        np_tprobs[i,] = leave_non_to_othernon
        # then row i at special state column is updated with leave_non_to_special.
        np_tprobs[i, special_state_idx] = leave_non_to_special
        # then row i at itself column is updated with stay_non.
        np_tprobs[i,i] = stay_non

    # Next, the matrix is updated at the special state row to reflect leave_special for all states.
    # from special state to any other state
    np_tprobs[special_state_idx,] = leave_special

    # Finally, the matrix is updated at the special state row and column to reflect stay_special for the special state self-transition.
    # from special state to itself
    np_tprobs[special_state_idx, special_state_idx] = stay_special

    ## normalize to 1
    for i in range(nstates):
        np_tprobs[i,] = np_tprobs[i,]/sum(np_tprobs[i,])
    
    return np_tprobs

def get_transition_probs_from_str(transprobstr):
    np_tprobs = np.array([[float(j) for j in e.split(',')] for e in transprobstr.split(';')])
    return np_tprobs

def get_transition_probs_exponential_decay(nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001, exp_decay_scale=1.0): ##, decay_rate=1):
    assert exp_decay_scale > 0
    leave_non = np.mean([leave_non_to_special, leave_non_to_othernon]) ## if --leave_other gave only one number then this is that number; if --leave_other gave two numbers then this is the average of those two numbers. 
    stay_special = 1 - (leave_special * (nstates-1))
    stay_non = 1 - leave_non ##leave_non_to_special - (leave_non_to_othernon * (nstates-2))
    np_tprobs = np.zeros([nstates,nstates], dtype=float)
    
    for i in range(nstates):
        for j in range(nstates):
            if i == special_state_idx:
                if j == special_state_idx:
                    np_tprobs[i,j] = stay_special
                else:
                    distance = abs(j - i) * exp_decay_scale
                    np_tprobs[i,j] = leave_special ** distance
            else:
                if j == i:
                    np_tprobs[i,j] = stay_non
                else:
                    distance = abs(j - i) * exp_decay_scale
                    # distance = (decay_rate ** (distance-1)) ## decay with distance from current state
                    np_tprobs[i,j] = leave_non_to_othernon ** distance

    ## normalize to 1
    for i in range(nstates):
        np_tprobs[i,] = np_tprobs[i,]/sum(np_tprobs[i,])
    return np_tprobs


def get_transition_probs(nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001, exp_decay=False, exp_decay_scale=1.0):
    
    if exp_decay:
        ## Here, leave probabilities decay with distance to the given state. So, for example, if you are in state 3, you are more likely to transition to state 4 than to state 7.
        np_tprobs = get_transition_probs_exponential_decay(nstates=nstates, special_state_idx=special_state_idx, leave_special=leave_special, leave_non_to_special=leave_non_to_special, leave_non_to_othernon=leave_non_to_othernon, exp_decay_scale=exp_decay_scale)
        ###################################################nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001
    else:
        np_tprobs = get_simple_transition_probs(nstates=nstates, special_state_idx=special_state_idx, leave_special=leave_special, leave_non_to_special=leave_non_to_special, leave_non_to_othernon=leave_non_to_othernon)
    
    return np_tprobs



    


def get_initial_probs(nstates=7, special_state_idx=0, special_state=0.997, other_states=None, inits=None):
    if inits is None:
        if other_states is None:
            other_states = (1.0-special_state)/float(nstates-1)
        inits = [other_states]*nstates
        inits[special_state_idx] = special_state
    np_iprobs = np.array(inits, dtype=float)
    np_iprobs = np_iprobs/np_iprobs.sum()
    return np_iprobs












## THESE ARE THE PROBABILITY MATRICES USED IN THE 7-STATE DNA PUFF FINDING MODEL

def emissions7():
    return get_emission_probs(mu=[1,2,4,8,16,32,64],sig=None)

def transitions7():
    return get_transition_probs(nstates=7, special_state_idx=0, leave_special=0.001, leave_non_to_special=0.001, leave_non_to_othernon=0.001)

def initial7():
    return get_initial_probs(nstates=7, special_state_idx=0, special_state=0.997)




##def ORI_bed(late):
##    ## BED with chr

def get_levels(states):
    levels = {}
    for chrom in states.keys():
        levels[chrom] = 2**(np.array(states[chrom])-1)
    return levels








### These help clean up the sub-program files
def help_get_emission_probs(mu, sigma=None, mu_scale=None, discrete=False):
    '''mu is a string w/ comma-separated state means
        sigma is a string with comma-sep stat stdevs
        RETURNS: e_prob matrix and nstates'''

    if discrete:
        e_mu = [[float(e) for e in E.split(',')] for E in mu.strip().strip('\\').split(';')]
    else:
        e_mu = [float(e) for e in mu.strip().strip('\\').split(',')]

    nstates = len(e_mu)

    if discrete:
        e_sig = None
    else:
        if sigma is not None:
            e_sig = [float(e) for e in sigma.strip().split(',')]
        elif mu_scale is not None: 
            e_sig = [e*mu_scale for e in e_mu]
        else:
            e_sig = [e**0.5 for e in e_mu]
        assert len(e_sig) == nstates

    np_eprobs = get_emission_probs(mu=e_mu, sig=e_sig, discrete=discrete)
    return np_eprobs, nstates


def help_get_transition_probs(leave_special_state, leave_other, special_state_idx, nstates, transprobs, exp_decay, exp_decay_scale):
    '''
        leave_special_state is probability of leaving special state (e.g. CN=1) - can make it same as others.
        leave_other is probability of leaving a non-special state
        special_state_idx is the 0-based idx of where to find special state params
        nstates is number of states in model
        transprobs = 
    '''
    if transprobs is not None:
        np_tprobs = get_transition_probs_from_str(transprobs)
    else:
        if leave_special_state > 1:
            leave_special_state = 1.0/leave_special_state
        if leave_other is None:
            leave_non_to_special = leave_special_state
            leave_non_to_other = leave_special_state
        else:
            leave_other = [float(e) for e in leave_other.split(',')]
            if len(leave_other) == 1:
                if leave_other[0] > 1:
                    lp = 1.0/leave_other[0]
                else:
                    lp = leave_other[0]
                leave_non_to_special = lp
                leave_non_to_other = lp
            elif len(leave_other) > 1:
                if leave_other[0] > 1:
                    leave_non_to_special = 1.0/leave_other[0]
                else:
                    leave_non_to_special = leave_other[0]
                if leave_other[1] > 1:
                    leave_non_to_other = 1.0/leave_other[1]
                else:
                    leave_non_to_other = leave_other[1]

        np_tprobs = get_transition_probs(nstates=nstates,
                                      special_state_idx=special_state_idx,
                                      leave_special=leave_special_state,
                                      leave_non_to_special= leave_non_to_special,
                                      leave_non_to_othernon=leave_non_to_other,
                                      exp_decay=exp_decay,
                                      exp_decay_scale=exp_decay_scale)
    return np_tprobs


def help_get_initial_probs(nstates, special_state_idx, init_special, initialprobs=None):
    if initialprobs is None:
        np_iprobs = get_initial_probs(nstates=nstates, special_state_idx=special_state_idx, special_state=init_special)
    else:
        inits = [float(e) for e in initialprobs.strip().split(',')]
        np_iprobs = get_initial_probs(inits=inits)
    assert len(np_iprobs) == nstates 
    return np_iprobs


def help_get_prob_matrices_from_params(mu, sigma, mu_scale, leave_special_state, leave_other, special_idx, init_special, initialprobs, transprobs, discrete=False, exp_decay=False, exp_decay_scale=1.0):
    np_eprobs, nstates = help_get_emission_probs(mu, sigma, mu_scale, discrete)
    np_tprobs = help_get_transition_probs(leave_special_state, leave_other, special_idx, nstates, transprobs, exp_decay=exp_decay, exp_decay_scale=exp_decay_scale)
    np_iprobs = help_get_initial_probs(nstates, special_idx, init_special, initialprobs)
    return np_eprobs, np_tprobs, np_iprobs






## KMEANS STUFF
def k_state_means(data, k):
    return kmeans(x=data, centers=k)

def get_k_mean_centers(km):
    ## km is kmeans object
    return list(km[1])

def get_k_mean_sigmas(km):
    ## km is kmeans object
    d = defaultdict(int)
    for k in km[0]:
        d[k-1] += 1 ## k-1 to put it in python indexing
    sigs = []
    for k in sorted(d.keys()):
        sig = (km[3][k]/d[k])**0.5
        sigs.append( sig )
    return sigs

def get_k_mean_cluster_lengths(km):
    d = defaultdict(list)
    last_c = km[0][0]
    L = 1
    for i in range(1,len(km[0])):
        this_c = km[0][i]
        if this_c == last_c:
            L+=1
        else:
            d[last_c].append(L)
            last_c = this_c
            L = 1
    d[this_c] = L
    return d ## all lengths found

def get_k_mean_cluster_mean_lengths(km_lengths):
    means = []
    for c in sorted(km_lengths.keys()):
        means.append( np.array(km_lengths[c]).mean() )
    return means

def get_k_mean_cluster_trans_probs(kmmeans):
    leave = 1/np.array(kmmeans, dtype=float)
    return list(1-(leave)), list(leave/(len(leave)-1))

def get_k_mean_cluster_initial_probs(km_lengths):
    sums = []
    for c in sorted(km_lengths.keys()):
        sums.append( np.array(km_lengths[c]).sum() )
    return list( np.array(sums, dtype=float)/np.array(sums,dtype=float).sum() )

## this method is not yet implemented in full work flow
def get_k_mean_cluster_trans_probs_matrix(km):
    t = np.ones([len(km[1]), len(km[1])], dtype=float) ## all start w/ pseudo count of 1
    for i in range(1,len(km[0])):
        c_from = km[0][i-1] - 1 #for python array indexing
        c_to = km[0][i] - 1 #for python array indexing
        t[c_from, c_to] += 1
    for i in range(len(km[1])):
        t[i,] = t[i,]/t[i,].sum()
    return t




def get_state_emissions_from_kmeans(data, k):
    km = k_state_means(data, k)
    mu = get_k_mean_centers(km)
    sig = get_k_mean_sigmas(km)
    km_lengths = get_k_mean_cluster_lengths(km)
    init_probs = get_k_mean_cluster_initial_probs(km_lengths)
    mean_lengths = get_k_mean_cluster_mean_lengths(km_lengths) #in order of clustnums like other2
    stay_probs, leave_probs = get_k_mean_cluster_trans_probs(mean_lengths)
    sorted_params = sorted( zip( mu, sig , init_probs, stay_probs, leave_probs) )
    ## outputs
    e_mu = [mu for mu,sig,I,S,L in sorted_params]
    e_sig = [sig for mu,sig,I,S,L in sorted_params]
    inits = [I for mu,sig,I,S,L in sorted_params]
    stay = [S for mu,sig,I,S,L in sorted_params]
    leave = [L for mu,sig,I,S,L in sorted_params]
    # returns
    return e_mu, e_sig, inits, stay, leave

def help_get_state_emissions_from_kmeans(data, k):
    raise NotImplementedError("kmeans functionality requires R/rpy2 which has been removed. Use alternative initialization.")
    # e_mu, e_sig, inits, stay, leave = get_state_emissions_from_kmeans(data, k)
    # np_eprobs = get_emission_probs(mu=e_mu, sig=e_sig)
    # rowsvec = []
    # for i in range(k):
    #     row = [leave[i]]*k
    #     row[i] = stay[i]
    #     rowsvec += row
    # np_tprobs = np.array(rowsvec).reshape((k,k))
    # np_iprobs = np.array(inits)
    # return np_eprobs, np_tprobs, np_iprobs
                                 
                            
## TODO: update above kmeans... to allow better transition estimates. e.g. can just make i,j matrix of all i-->j
























## FUNCTIONS NOT BEING USED
def window_medians(x, flanksize=2, includeflanks=True):
    ''' Given vector of FEs and flank size, find median FE in window centered over each bin by entending flanksize in each direction.
    TODO: For edge cases, take median in windows of flanksize*2 by adding/subtracting on each side as nec...'''
    meds = []
    lenx = len(x)
    start = flanksize
    end = lenx-flanksize
    Pos = start
    ans = np.median( x[Pos-flanksize:Pos+1+flanksize] )
    if includeflanks:
        # Append ans for all flank positions before Start
        meds += [ans]*flanksize
    meds.append(ans)
    #Iter from start+1
    for Pos in range(start+1, end):
        meds.append( np.median( x[Pos-flanksize:Pos+1+flanksize] ) )
    if includeflanks:
        # Append ans for all flank positions before Start
        meds += [meds[-1]]*flanksize
    return np.array(meds)
    

def window_sums(x, flanksize=2, returnmeans=False, includeflanks=True):
    ''' Given vector and flank size, get sum or mean of bin and flanks to each side.
        Note that this is sliding window mean-smoothing with uniform weighting.
        Can use ksmooth or loess in R to weight closer bins higher.'''
    sums = []
    lenx = len(x)
    start = flanksize
    end = lenx-flanksize 
    windowsize = float(flanksize*2 + 1) ## 2 flanks + Pos it is centered on
    Pos = start
    #First window
    ans = np.sum( x[Pos-flanksize:Pos+flanksize+1] )
    if includeflanks:
        # Append ans for all flank positions before Start
        sums += [ans]*flanksize
    # Append ans for Start
    sums.append( ans ) 

    #Iterate
    for Pos in range(start+1, end):
        ## Subtract first element of last window
        ans -= x[Pos-1-flanksize]
        ## Add last element of current window
        ans += x[Pos+flanksize]
        ## Append ans to sums
        sums.append( ans )

    if includeflanks:
        # Append last answer for all flank positions after End
        sums += [sums[-1]] * flanksize
        
    ## convert sums to np
    sums = np.array(sums)
    ## If means desired
    if returnmeans:
        return sums/windowsize
    return sums


def ksmooth_counts(self, bw=10000):
    raise NotImplementedError("ksmooth_counts requires R/rpy2 which has been removed. Use scipy alternatives.")
    # for chrom in self.chromosomes:
    #     x = self.start[chrom]
    #     y = self.count[chrom]
    #     k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
    #     self.count[chrom] = np.array(k[1])

## HMM followed by window_modes of states could help...
