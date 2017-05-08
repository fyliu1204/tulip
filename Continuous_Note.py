
def discretBisimulation(

    part, sys, N=10, min_cell_volume=0.1,

    closed_loop=True,

    max_num_poly=5, use_all_horizon=False,

    trans_length=1, remove_trans=False,

    abs_tol=1e-7

):

"""

#for discrete system, look at synth.synthesize for example (sys)

    @param part: PropPreservingPartition

    @param sys: discrete system

    @param N: horizon length

    @param min_cell_volume: the minimum volume of cells in the resulting

        partition.

    @param closed_loop: boolean indicating whether the `closed loop`

        algorithm should be used. default True.

    @param max_num_poly: maximum number of polytopes in a region to use in

        reachability analysis.

    @param use_all_horizon: in closed loop algorithm: if we should look

        for reachability also in less than N steps.

    @param trans_length: the number of polytopes allowed to cross in a

        transition.  a value of 1 checks transitions

        only between neighbors, a value of 2 checks

        neighbors of neighbors and so on.

    @param remove_trans: if True, remove found transitions between

        non-neighbors.

    @param abs_tol: maximum volume for an "empty" polytope

    """





    min_cell_volume = (min_cell_volume /np.finfo(np.double).eps

        *np.finfo(np.double).eps)



    # step 2. initialize transition matrix TR:

    n = len(part)

    transitions = np.zeros([n, n], dtype = int)



    #step 4. neighbor matrix NB

    IJ = part.adj.copy()

    IJ = IJ.todense()

    IJ = np.array(IJ)



    # step 5. def reachable_within(trans_length, adj_k, adj):

    IJ = reachable_within(trans_length, IJ,

                          np.array(part.adj.todense()) )



    # Initialize output, adj:copy of NB

    sol = deepcopy(part.regions) # step 3. set S

    adj = part.adj.copy()

    adj = adj.todense()

    adj = np.array(adj)





    progress = list()

    ss = sys



    # Step 6. Do the abstraction

    while np.sum(IJ) > 0:



        #step7. Pick (i; j) such that IJ( j; i) = 1, then set IJ( j; i) = 0.

        ind = np.nonzero(IJ)

        i = ind[1][0]

        j = ind[0][0]

        IJ[j, i] = 0



        # si, sj belongs to S

        si = sol[i]

        sj = sol[j]





        si_tmp = deepcopy(si)

        sj_tmp = deepcopy(sj)





        # step 8. solve_feasible

        S0 = solve_feasible(

            si, sj, ss, N, closed_loop,

            use_all_horizon, trans_set, max_num_poly

        )





        #intersect

        isect = si.intersect(S0)

        vol1 = isect.volume

        risect, xi = pc.cheby_ball(isect)



        #diff

        diff = si.diff(S0)

        vol2 = diff.volume

        rdiff, xd = pc.cheby_ball(diff)



        if vol1 <= min_cell_volume:

            logger.warning('\t too small: si \cap Pre(sj), ' +

                           'so discard intersection')

        if vol1 <= min_cell_volume and isect:

            logger.warning('\t discarded non-empty intersection: ' +

                           'consider reducing min_cell_volume')

        if vol2 <= min_cell_volume:

            logger.warning('\t too small: si \ Pre(sj), so not reached it')



        # We don't want our partitions to be smaller than the disturbance set

        # Could be a problem since cheby radius is calculated for smallest

        # convex polytope, so if we have a region we might throw away a good

        # cell.



        # step 9. if (vol(Vi \ S0) > min_vol) AND (vol(VinS0) > min_vol) then

        if (vol1 > min_cell_volume) and (risect > rd) and \

           (vol2 > min_cell_volume) and (rdiff > rd):



            # Make sure new areas are Regions and add proposition lists

            if len(isect) == 0:

                isect = pc.Region([isect], si.props)

            else:

                isect.props = si.props.copy()



            if len(diff) == 0:

                diff = pc.Region([diff], si.props)

            else:

                diff.props = si.props.copy()



            # step 10. replace si by intersection (single state)

            """

            def separate(reg1, abs_tol=ABS_TOL): Divide a region into several

            regions such that they are all connected.

            @type reg1: L{Region}

            @param abs_tol: Absolute tolerance

            @return: List [] of connected Region

            """

            isect_list = pc.separate(isect)

            sol[i] = isect_list[0]



            # cut difference into connected pieces ???

            difflist = pc.separate(diff)

            difflist += isect_list[1:]

            n_isect = len(isect_list) -1

            num_new = len(difflist)



            # step 11

            for region in difflist:

                sol.append(region)

            n_cells = len(sol)

            new_idx = xrange(n_cells-1, n_cells-num_new-1, -1)



            # step 12. Add row and column of zeros to TR

            transitions = np.pad(transitions, (0,num_new), 'constant')

            transitions[i, :] = np.zeros(n_cells)

            for r in new_idx: #why do we need this?

                transitions[i, r] = 0

                transitions[j, r] = 0





            # Step 15-17

            if i != j:

                transitions[j, i] = 1





            # step 12-13. Update NB(end; :), NB(:;end), NB(i; :), NB(:; i)



            old_adj = np.nonzero(adj[i, :])[0] #index of i-th row which is nonzero

            adj[i, :] = np.zeros([n_cells - num_new])

            adj[:, i] = np.zeros([n_cells - num_new])

            adj[i, i] = 1



            adj = np.pad(adj, (0, num_new), 'constant')



            for r in new_idx:

                adj[i, r] = 1

                adj[r, i] = 1

                adj[r, r] = 1



            for r in new_idx:

                for k in new_idx:

                    if r is k:

                        continue

                    if pc.is_adjacent(sol[r], sol[k]):

                        adj[r, k] = 1

                        adj[k, r] = 1





            for k in np.setdiff1d(old_adj, [i,n_cells-1]):

                if pc.is_adjacent(sol[i], sol[k]):

                    adj[i, k] = 1

                    adj[k, i] = 1

                elif remove_trans and (trans_length == 1):

                    transitions[i, k] = 0

                    transitions[k, i] = 0



                for r in new_idx:

                    if pc.is_adjacent(sol[r], sol[k]):

                        adj[r, k] = 1

                        adj[k, r] = 1

                    elif remove_trans and (trans_length == 1):

                        transitions[r, k] = 0

                        transitions[k, r] = 0



            # step 18-22. Update IJ matrix

            IJ = np.pad(IJ, (0,num_new), 'constant')

            adj_k = reachable_within(trans_length, adj, adj)

            sym_adj_change(IJ, adj_k, transitions, i)



            for r in new_idx:

                sym_adj_change(IJ, adj_k, transitions, r)



        # step 23-24:

        elif vol2 < abs_tol:

            transitions[j,i] = 1



    output = {'partition': sol, 'transition': transitions}

    return output





def sym_adj_change(IJ, adj_k, transitions, i):

    horizontal = adj_k[i, :] -transitions[i, :] > 0

    vertical = adj_k[:, i] -transitions[:, i] > 0



    IJ[i, :] = horizontal.astype(int)

    IJ[:, i] = vertical.astype(int)



def reachable_within(trans_length, adj_k, adj):

    """Find cells reachable within trans_length hops.

    """

    if trans_length <= 1:

        return adj_k



    k = 1

    while k < trans_length:

        adj_k = np.dot(adj_k, adj)

        k += 1

    adj_k = (adj_k > 0).astype(int)



    return adj_k