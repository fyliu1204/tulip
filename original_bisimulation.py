from scipy.sparse import lil_matrix
from numpy.random import rand
import networkx as nx
import tulip.transys as trs
import warnings
'''
Function: coarsest possibel proposition preserving partition
Input:  A finite transition system
Output: A list. Each entry of the list is a subset of states set S. Elements in
        the same entry have the same atomic proposition. Elements in different 
        entries have different atomic propositions.
'''
def coarsest_possible_ppp(fts):
    d1 = dict()
    for states in fts.states():
        ap = fts.states.find(states)
        (key, value) = ap[0]
        ap_list = repr(value['ap'])
        if ap_list not in d1:
            d1[ap_list] = []
        d1[ap_list].append(states)
    c_list = list()
    for key in d1:
        c_list.append(d1[key])
    return c_list
'''
Function: Original Bisimulation Algorithm
It's the implementation of Algorithm 1 in
http://web.eecs.umich.edu/~necmiye/pubs/WagenmakerO_allerton16.pdf
Input:  A finite transition system
Output: Another finite transition system. It is the bisimulation system of the
        original finite transition system
'''
def bisimulation_algorithm(fts):
    c_list = coarsest_possible_ppp(fts)
    flag = True
    while flag:
        flag = False
        for si in c_list:
            if flag == True:
                break
            for sj in c_list:
                pre_sj = list(fts.states.pre(si))
                tmp = [val for val in si if val in pre_sj]
                if len(tmp) != 0 and tmp != si:
                    flag = True
                    c_list.remove(si)
                    c_list.append(tmp)
                    c_list.append(list(set(si).difference(set(pre_sj))))
                    break
    new_fts = trs.FTS()
    #still working on it
    return new_fts