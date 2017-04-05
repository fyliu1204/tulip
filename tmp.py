Step 2:
- To find the ap of a single state C{'s0'}:
              >>> a = ts.states.find(['s0'] )
              >>> (s0_, label) = a[0]
              >>> print(label)
              {'ap': set(['p'])}
- To find all states with a specific ap C{{'p'}}:
              >>> ts.states.add('s1', ap={'p'})
              >>> b = ts.states.find(with_attr_dict={'ap':{'p'} } )
              >>> states = [state for (state, label_) in b]
              >>> print(set(states) )
              {'s0', 's1'}

              d1 = dict{}
              for states in fts.states():
                  ap = fts.states.find(states)
                  (key, value) = ap[0]
                  if value['ap'] not in d1:
                      d1[value['ap']] = []
                  d1[value['ap']].append(states)

labeled_graphs.py States find


Step 3:
"""
https://github.com/tulip-control/tulip-control/blob/053101e8e743c87585daf33ddb479e7bfa0dcfd6/tulip/transys/labeled_graphs.py
line 232:
def pre(self, states):
        """Return direct predecessors (1-hop) of given state.
        See Also
        ========
          - L{post}
          - Def. 2.3, p.23 U{[BK08]
            <https://tulip-control.sourceforge.io/doc/bibliography.html#bk08>}
        @rtype: set
        """
        states = self._single_state2singleton(states)
        predecessors = set()
        for state in states:
            predecessors |= set(self.graph.predecessors(state))
        return predecessors
e.g. pre_all = fts.states.pre(['select','soda'])
input_type: list return type:set
"""


Questions:
1. input type?
fts = trs.FTS()
2. Algorithm 1: s1, s2 type? subset of states?
Type of S?
