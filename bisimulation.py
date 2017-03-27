def bisimulation(sys):

    # @param sys: transition system object


    # Assume S as coarsest possible proposition preserving partition
    # coarsest possible proposition preserving partition algorithm not finished
    n = len(S)
    transitions = np.zeros([n, n], dtype = int)
    IJ = part.adj.copy()
    IJ = IJ.todense()
    IJ = np.array(IJ)
