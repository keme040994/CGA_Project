# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: April, 2016.

# LIBRARIES
import numpy as np
import scipy.stats.stats as sss


# ==============
# DATA FUNCTIONS
# ==============

# FUNCTION: data_switcher
def data_switcher(a_switch_log, a_switch_zscore, rep):
    """
    An auxiliary function that receives some flags and the replications to determine if the log() and zscore() must be
    applied.

    Args:
        a_switch_log : BOOLEAN
            A flag value for the log() function
        a_switch_zscore : BOOLEAN
            A flag value for the zscore() function
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result

    Returns:
         LIST[rep1, rep2, rep3, ...]
              The new transformed replications
    """
    if a_switch_log:
        rep = log_function(rep)
    if a_switch_zscore:
        rep = zscore_function(rep)
    return rep


# FUNCTION: log
def log_function(rep):
    """
    This function applies the log() transform to every value in the replications.

    Args:
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result

    Returns:
         LIST[rep1, rep2, rep3, ...]
              The new transformed replications
    """
    rep = map(lambda x: np.matrix([[np.log(item) for item in row] for row in x.tolist()]), rep)
    return list(rep)


# FUNCTION: zcore
def zscore_function(rep):
    """
    This function applies the zscore() transform to every value in the replications.

    Args:
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result

    Returns:
         LIST[rep1, rep2, rep3, ...]
              The new transformed replications
    """
    rep = map(lambda x: np.asmatrix(np.transpose(np.array([sss.zscore(item) for item in np.transpose(np.asarray(x))]))),
              rep)
    return list(rep)
