from tulip.transys import labeled_graphs
import polytope as pc

def bisimulation(sys):

    # @param sys: transition system object

    """
    Assume part as coarsest possible proposition preserving partition
    coarsest possible proposition preserving partition algorithm not finished
    part: class PropPreservingPartition(pc.MetricPartition), in prop2partition
    pre: exists in labeled_graphs, rtype: set()
    region: list of polytopes
    """
    n = len(part)
    transitions = np.zeros([n, n], dtype = int)


    S = deepcopy(part.regions) #type: list of region

    data = checkCondition(S)
    while (data['exists'] == True):
        i = data['i']
        j = data['j']
        si = S[i] #type: region
        sj = S[j]
        isect = data['isect'] # si isect pre(sj)
        diff = data['diff'] # si diff pre(sj)
        del S[i] # S/si
        S = S.append(isect)
        S = S.append(diff)
        data = checkCondition(S)




def checkCondition(S):
    n = length(S)
    exists = False
    for i in range(0, n):
        si = S[i] # type: region
        for j in range(0, n):
            sj = S[j] # type: region
            # can i equal to j?
            pred = pre(j) # set of numbers
            """
            not sure if this union and diff is correct:
            intersect: line 603 in polytope.py, rtype: region
            """
            r1 = Region()
            r2 = Region()
            for k in pred:
                tmp = S[k] # type: region
                isect = si.intersect(tmp) # rtype: region
                diff = si.diff(tmp) # rtype: region
                r1 = union(r1, isect, check_convex = True)
                r2 = union(r2, diff, check_convex = True)
            vol1 = r1.volume
            vol2 = r2.volume

            if (vol1 > 0) and (vol2 > 0):
                exists = True
                data = {
                    'exists': exists,
                    'isect': r1,
                    'diff': r2,
                    'i': i,
                    'j': j
                }
                return data

    data = {
        'exists': exists,
        'isect': None,
        'diff': None,
        'i': None,
        'j': None
    }
    return data
