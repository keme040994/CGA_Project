# Created by: Dr. David John
# July 2015:
#      Original version
# February 2016:
#      Significant improvement to commenting
# April 2016
#      - Fixed a division problem in computing the
#      average of the parents value (1/_r) with (1/float(_r))
#      - Fixed division by 2 in answer3 to be division by 2.0

# LIBRARIES
import functools
import math
import networkx as nx
import numpy as np


# pnProb(G,rep)
# compute the cotemporal marginal likelihood of DAG G given replicate
# set rep using the Patton-Norris marginal likelihood formula

# FUNCTION: pnml_cotemporal
def pnml_cotemporal(di_graph, rep):
    """
    Given a DAG G and a set of replicates, rep, compute the cotemporal
    marginal likelihood
    that G describes the data, pnML(G,rep).  For details on this likelihood
    see pages 27-33 of Kris Patton's master's thesis, Department of
    Mathematics, 2012.  Terms g, v0, sigma are all equal to 1.0.
    """
    # number of time points and number of replicates
    _t = np.shape(rep)[1]
    _r = len(rep)

    # partition the nodes of G into two classes, those without (wo) parents and
    # those with at least one parent
    _VwoP = [_v for _v in nx.nodes(di_graph) if nx.ancestors(di_graph, _v) == set([])]
    _VP = [_v for _v in nx.nodes(di_graph) if _v not in _VwoP]

    # compute (and show) likelihood contribution for nodes with out parents
    _MLwithOutParents = math.pow(2.0*math.pi*math.e,
                                 -0.5*_r*_t*len(_VwoP))/math.e

    # compute (and show) likelihood contribution for nodes with at least one parent

    # dict that will be used to precompute a number of values and matrices that are used
    # in the later computation
    _nodeInfo = {}
    for _child in _VP:
        # get all time data for child (cotemporal paradigm)
        _child_data = [_zzz[:, _child] for _zzz in rep]

        # get parents for child and all parent's time data (cotemporal paradigm)
        _Parents = [_parent for _parent in nx.ancestors(di_graph, _child)
                    if di_graph.has_edge(_parent, _child)]
        _parent_data = [_zzz[:, _Parents] for _zzz in rep]

        # compute average of the parents data over the replicates t X Kc
        if _r == 1:
            __avg_parent_data = _parent_data[0]
        else:
            __avg_parent_data = (1/float(_r))*functools.reduce(lambda x, y: x+y, _parent_data)

        # There are a number of values that used in numerous places
        # in the calculation of the likelihood.  Some of these are
        # fairly computationally expensive, hence  the strategy is to
        # compute these once, store them in a data dictionary, and then
        # retrieve these precomputed values as needed in the final
        # calculation.  For each of the "r" replicates there are three
        # such values that will be precomputed, and there is one final
        # value as well. Hence "3*r+1" values will be stored in the
        # data store "nodeInfo".

        # For specific details see the closed form expression of
        # f(Q | DAG) on page 33 of Patton's master's thesis.

        # r pieces of data for child
        _holder = [np.dot(np.transpose(_v), _v) for _v in _child_data]

        # r pieces of data associated with parents of child
        _holder.extend([np.dot(np.transpose(_v), _v)
                        for _v in _parent_data])

        # r pieces of data associated with child and parents
        _holder.extend([np.dot(np.transpose(_v), _w)
                        for _v, _w in zip(_parent_data, _child_data)])

        # averages over all replicates of parents
        _holder.append(np.dot(np.transpose(__avg_parent_data), __avg_parent_data))

        # put in dictionary
        _nodeInfo[_child] = _holder

    # Compute the Patton-Norris marginal likelihood using
    # all the precomputed pieces.  For the details of this
    # computations go see Kris Patton's thesis, page 33.
    # This computation will use the precomputed values
    # stored in nodeInfo.

    # initialize floating point working variables
    _answer1 = _answer2 = _answer3 = 1.0

    for _val in _nodeInfo.values():
        # product across all nodeInfo's of inner products
        # of avg. of parents raised to "r/2" power
        _answer1 *= math.pow(np.linalg.det(_val[3*_r]),0.5*_r)

        # product across all nodeInfo's of square root of
        # determinant of the sum of inner product of parents and
        # inner product of avg. of parents
        for _index in range(_r, 2*_r):
            _answer2 *= math.sqrt(
                np.linalg.det(_val[_index]+_val[3*_r]))

        # xxxtemp is the sum across all replicates, initially a "1"
        _xxxtemp = 1.0
        for _index in range(0, _r):
            # Values used on _xxxtemp
            # a) inner product of child
            # b) transpose(A)
            # c) inverse(inner product of parent - inner product avg parent
            # d) A
            _xxxtemp += \
                    _val[_index] - \
                    np.transpose(_val[2*_r+_index]) * \
                    np.transpose(np.linalg.inv(
                        (_val[_r+_index]+_val[3*_r]))) * \
                    _val[2*_r+_index]

        # xxxtemp raised to (-r*t+1)/2 power
        _answer3 *= math.pow(_xxxtemp, -(_r*_t+1)/2.0)

    # putting all the three major pieces together
    _MLwithParents = _answer3*_answer1/_answer2

    # coef of the likelihood term
    _MLcoef = math.pow(math.pi, -_r*_t*len(_VP)/2.0) * math.pow(math.gamma((_r*_t+1)/2.0), len(_VP))

    return _MLcoef*_MLwithParents*_MLwithOutParents


# FUNCTION: pnml_next_step_one
def pnml_next_step_one(di_graph, rep):
    # number of time points and number of replicates
    _t = np.shape(rep)[1]
    _r = len(rep)

    # partition the nodes of G into two classes, those without (wo) parents and
    # those with at least one parent
    _VwoP = [_v for _v in nx.nodes(di_graph) if nx.ancestors(di_graph, _v) == set([])]
    _VP = [_v for _v in nx.nodes(di_graph) if _v not in _VwoP]

    # compute (and show) likelihood contribution for nodes with out parents
    _MLwithOutParents = math.pow(2.0*math.pi*math.e,
                                 -0.5*_r*(_t-1)*len(_VwoP))/math.e

    # compute (and show) likelihood contribution for nodes with at least one parent

    # dict that will be used to precompute a number of values and matrices that are used
    # in the later computation
    _nodeInfo = {}
    for _child in _VP:
        # get all time data for child ('next step one' paradigm)
        _child_data = [_zzz[1:(_t-1), _child] for _zzz in rep]

        # get parents for child and all parent's time data ('next step one' paradigm)
        _Parents = [_parent for _parent in nx.ancestors(di_graph, _child)
                    if di_graph.has_edge(_parent, _child)]
        _parent_data = [_zzz[0:(_t-2), _Parents] for _zzz in rep]

        # compute average of the parents data over the replicates (t-1) X Kc
        if _r == 1:
            __avg_parent_data = _parent_data[0]
        else:
            __avg_parent_data = (1/float(_r))*functools.reduce(lambda x, y: x+y, _parent_data)

        # There are a number of values that used in numerous places
        # in the calculation of the likelihood.  Some of these are
        # fairly computationally expensive, hence  the strategy is to
        # compute these once, store them in a data dictionary, and then
        # retrieve these precomputed values as needed in the final
        # calculation.  For each of the "r" replicates there are three
        # such values that will be precomputed, and there is one final
        # value as well. Hence "3*r+1" values will be stored in the
        # data store "nodeInfo".

        # For specific details see the closed form expression of
        # f(Q | DAG) on page 33 of Patton's master's thesis.

        # r pieces of data for child
        _holder = [np.dot(np.transpose(_v), _v) for _v in _child_data]

        # r pieces of data associated with parents of child
        _holder.extend([np.dot(np.transpose(_v), _v)
                        for _v in _parent_data])

        # r pieces of data associated with child and parents
        _holder.extend([np.dot(np.transpose(_v), _w)
                        for _v, _w in zip(_parent_data, _child_data)])

        # averages over all replicates of parents
        _holder.append(np.dot(np.transpose(__avg_parent_data), __avg_parent_data))

        # put in dictionary
        _nodeInfo[_child] = _holder

    # Compute the Patton-Norris marginal likelihood using
    # all the precomputed pieces.  For the details of this
    # computations go see Kris Patton's thesis, page 33.
    # This computation will use the precomputed values
    # stored in nodeInfo.

    # initialize floating point working variables
    _answer1 = _answer2 = _answer3 = 1.0

    for _val in _nodeInfo.values():
        # product across all nodeInfo's of inner products
        # of avg. of parents raised to "r/2" power
        _answer1 *= math.pow(np.linalg.det(_val[3*_r]),0.5*_r)

        # product across all nodeInfo's of square root of
        # determinant of the sum of inner product of parents and
        # inner product of avg. of parents
        for _index in range(_r, 2*_r):
            _answer2 *= math.sqrt(
                np.linalg.det(_val[_index]+_val[3*_r]))

        # xxxtemp is the sum across all replicates, initially a "1"
        _xxxtemp = 1.0
        for _index in range(0, _r):
            # Values used on _xxxtemp
            # a) inner product of child
            # b) transpose(A)
            # c) inverse(inner product of parent - inner product avg parent
            # d) A
            _xxxtemp += \
                    _val[_index] - \
                    np.transpose(_val[2*_r+_index]) * \
                    np.transpose(np.linalg.inv(
                        (_val[_r+_index]+_val[3*_r]))) * \
                    _val[2*_r+_index]

        # xxxtemp raised to (-r*t+1)/2 power
        _answer3 *= math.pow(_xxxtemp, -(_r*(_t-1)+1)/2.0)

    # putting all the three major pieces together
    _MLwithParents = _answer3*_answer1/_answer2

    # coef of the likelihood term
    _MLcoef = math.pow(math.pi, -_r*(_t-1)*len(_VP)/2.0) * math.pow(math.gamma((_r*(_t-1)+1)/2.0), len(_VP))

    return _MLcoef*_MLwithParents*_MLwithOutParents


# FUNCTION: pnml_next_step_one_two
def pnml_next_step_one_two(di_graph, rep):
    # number of time points and number of replicates
    _t = np.shape(rep)[1]
    _r = len(rep)

    # partition the nodes of G into two classes, those without (wo) parents and
    # those with at least one parent
    _VwoP = [_v for _v in nx.nodes(di_graph) if nx.ancestors(di_graph, _v) == set([])]
    _VP = [_v for _v in nx.nodes(di_graph) if _v not in _VwoP]

    # compute (and show) likelihood contribution for nodes with out parents
    _MLwithOutParents = math.pow(2.0*math.pi*math.e,
                                 -0.5*_r*(_t-2)*len(_VwoP))/math.e

    # compute (and show) likelihood contribution for nodes with at least one parent

    # dict that will be used to precompute a number of values and matrices that are used
    # in the later computation
    _nodeInfo = {}
    for _child in _VP:
        # get all time data for child ('next step one-two' paradigm)
        _child_data = [_zzz[2:(_t-1), _child] for _zzz in rep]

        # get parents for child and all parent's time data ('next step one-two' paradigm)
        _Parents = [_parent for _parent in nx.ancestors(di_graph, _child)
                    if di_graph.has_edge(_parent, _child)]
        _parent_data = [_zzz[0:(_t-3), _Parents] for _zzz in rep]

        # compute average of the parents data over the replicates (t-2) X Kc
        if _r == 1:
            __avg_parent_data = _parent_data[0]
        else:
            __avg_parent_data = (1/float(_r))*functools.reduce(lambda x, y: x+y, _parent_data)

        # There are a number of values that used in numerous places
        # in the calculation of the likelihood.  Some of these are
        # fairly computationally expensive, hence  the strategy is to
        # compute these once, store them in a data dictionary, and then
        # retrieve these precomputed values as needed in the final
        # calculation.  For each of the "r" replicates there are three
        # such values that will be precomputed, and there is one final
        # value as well. Hence "3*r+1" values will be stored in the
        # data store "nodeInfo".

        # For specific details see the closed form expression of
        # f(Q | DAG) on page 33 of Patton's master's thesis.

        # r pieces of data for child
        _holder = [np.dot(np.transpose(_v), _v) for _v in _child_data]

        # r pieces of data associated with parents of child
        _holder.extend([np.dot(np.transpose(_v), _v)
                        for _v in _parent_data])

        # r pieces of data associated with child and parents
        _holder.extend([np.dot(np.transpose(_v), _w)
                        for _v, _w in zip(_parent_data, _child_data)])

        # averages over all replicates of parents
        _holder.append(np.dot(np.transpose(__avg_parent_data), __avg_parent_data))

        # put in dictionary
        _nodeInfo[_child] = _holder

    # Compute the Patton-Norris marginal likelihood using
    # all the precomputed pieces.  For the details of this
    # computations go see Kris Patton's thesis, page 33.
    # This computation will use the precomputed values
    # stored in nodeInfo.

    # initialize floating point working variables
    _answer1 = _answer2 = _answer3 = 1.0

    for _val in _nodeInfo.values():
        # product across all nodeInfo's of inner products
        # of avg. of parents raised to "r/2" power
        _answer1 *= math.pow(np.linalg.det(_val[3*_r]),0.5*_r)

        # product across all nodeInfo's of square root of
        # determinant of the sum of inner product of parents and
        # inner product of avg. of parents
        for _index in range(_r, 2*_r):
            _answer2 *= math.sqrt(
                np.linalg.det(_val[_index]+_val[3*_r]))

        # xxxtemp is the sum across all replicates, initially a "1"
        _xxxtemp = 1.0
        for _index in range(0, _r):
            # Values used on _xxxtemp
            # a) inner product of child
            # b) transpose(A)
            # c) inverse(inner product of parent - inner product avg parent
            # d) A
            _xxxtemp += \
                    _val[_index] - \
                    np.transpose(_val[2*_r+_index]) * \
                    np.transpose(np.linalg.inv(
                        (_val[_r+_index]+_val[3*_r]))) * \
                    _val[2*_r+_index]

        # xxxtemp raised to (-r*t+1)/2 power
        _answer3 *= math.pow(_xxxtemp, -(_r*(_t-2)+1)/2.0)

    # putting all the three major pieces together
    _MLwithParents = _answer3*_answer1/_answer2

    # coef of the likelihood term
    _MLcoef = math.pow(math.pi, -_r*(_t-2)*len(_VP)/2.0) * math.pow(math.gamma((_r*(_t-2)+1)/2.0), len(_VP))

    return _MLcoef*_MLwithParents*_MLwithOutParents
