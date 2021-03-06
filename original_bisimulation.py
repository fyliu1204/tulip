'''
这是我为Python for The Temporal Logic Planning (TuLiP) toolbox加入的几个函数和几个简单的testcase。
在编写函数中，大量使用了之前已有的api和数据结构， 并将新函数很好的嵌入到该toolbox中。
'''
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
Function: List of states transfer to a fts
Input:  A list. Each entry of the list is a subset of states set S. Elements in
        the same entry have the same atomic proposition. Elements in different 
        entries have different atomic propositions.
        A finite transition system
Output: Another finite transition system. It is the bisimulation system of the
        original finite transition system
Comment: We ignored the environment actions and system actions
'''
def list_to_fts(c_list, fts):
    print c_list
    new_fts = trs.FTS()
    new_fts.name = fts.name
    for state in c_list:
        a = fts.states.find([state[0]])
        (s0_, label) = a[0]
        new_fts.states.add(''.join(state))
        for ap in label['ap']:
            new_fts.atomic_propositions.add(ap)
        new_fts.states[''.join(state)]['ap'] = label['ap']
    for sys_actions in fts.sys_actions:
        new_fts.sys_actions.add(sys_actions)
    for env_actions in fts.env_actions:
        new_fts.env_actions.add(env_actions)
    for ini_state in fts.states.initial:
        if ini_state in new_fts.states:
            new_fts.states.initial.add(ini_state)
        else:
            for state in new_fts.states:
                if ini_state in state:
                    new_fts.states.initial.add(state)
                    break
    for tran in fts.transitions(data=True):
        s1 = ''
        s2 = ''
        if tran[0] in new_fts.states:
            s1 = tran[0]
        else:
            for state in new_fts.states:
                if tran[0] in state:
                    s1 = state
                    break
        if tran[1] in new_fts.states:
            s2 = tran[1]
        else:
            for state in new_fts.states:
                if tran[1] in state:
                    s2 = state
                    break
        flag = False
        for tmp in new_fts.transitions():
            if s1 == tmp[0] and s2 == tmp[1]:
                flag = True
                break
        if flag:
            continue
        sys_act = ''
        env_act = ''
        if 'sys_actions' in tran[2]:
            sys_act = tran[2]['sys_actions']
        if 'env_actions' in tran[2]:
            env_act = tran[2]['env_actions']
        new_fts.transitions.add(
            s1, s2
        )
    print new_fts
    return new_fts
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
    new_fts = list_to_fts(c_list, fts)
    return new_fts

def dual_simulation_algorithm(fts):
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
                if len(tmp) != 0 and tmp not in c_list:
                    flag = True
                    c_list.append(tmp)
                    break
    new_fts = list_to_fts(c_list, fts)
    return new_fts

ofts = trs.FiniteTransitionSystem()

ofts.states.add_from(['s1', 's2', 's3', 's4'] )
ofts.states.initial.add('s1')

ofts.atomic_propositions |= ['p1','p2']
apset = set()
apset.add('p1')
apset.add('p2')
ofts.states['s1']['ap'] = apset
ofts.states.add('s2', ap={'p1'} )
ofts.states.add('s3', ap={'p2'})
ofts.states.add('s4', ap={'p1', 'p2'})

ofts.transitions.add('s1', 's4') # unlabeled

ofts.sys_actions.add('try')
ofts.sys_actions.add_from({'start', 'stop'} )
ofts.env_actions.add_from({'block', 'wait'} )

# remove unlabeled edge s1->s2
#ofts.transitions.remove('s1', 's2')

ofts.transitions.add(
    's2', 's1',
    sys_actions='try', env_actions='block'
)
ofts.transitions.add(
    's3', 's2',
    sys_actions='start', env_actions='wait'
)
ofts.transitions.add(
    's4', 's3',
    sys_actions='stop', env_actions='block'
)
print ofts
fts = dual_simulation_algorithm(ofts)
