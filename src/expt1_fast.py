import cplex
import sys
import os
import networkx as nx
import numpy as np
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from cplex.exceptions import CplexError
from time import perf_counter
import matplotlib.pyplot as plt
import numpy
import time
import chrvrp_io

LPlowerbound = -1_000_000  # was `-10 ^ 6`: `^` is XOR in Python, not power, so this evaluated to -16
### RUN PARAMETERS ###
# --- perf/debug switches (added during cleanup) ---
DEBUG_CUTS = 0    # 1 = re-enable the per-arc diagnostic prints inside the callback
SPARSE_CUTS = 0   # 1 = drop zero coefficients from Benders cuts (see notes)
multicut = 0
tight2 = 0
table = 0
# Experiment switches, read from the environment so batch runs need no edits:
#   CHRVRP_LB="1234"          lower bounds to add (any subset of 1-4; "" = none)
#   CHRVRP_BETATIGHT=1        beta-tightening of tailored Benders cuts
#   CHRVRP_FORCED_VEHICLES=1  C4 "= 1" (every vehicle visits a supplier);
#                             default is the paper's (da5) "<= 1"
_lbs = os.environ.get('CHRVRP_LB', '')
lowerb1 = int('1' in _lbs)
lowerb2 = int('2' in _lbs)
lowerb3 = int('3' in _lbs)
lowerb4 = int('4' in _lbs)
betatight = int(os.environ.get('CHRVRP_BETATIGHT', 0))
FORCED_VEHICLES = int(os.environ.get('CHRVRP_FORCED_VEHICLES', 0))
gap = 0
bound =-1

#filename1 = "instance10S2V3Y4PxLxn29.tsp"
#capacitym = 6
#runcode = 1

# --- named run mode ---
# Each name fully determines both the solution method (runcode, used by
# createILP()/solve_and_report() below) and which link constraint defines z
# (MODEL -> zxlink/zflink), so there is no separate "which model" argument
# needed -- picking the mode also picks the model.
#   zxdefault  -> D_z-x, compact formulation, default CPLEX solver     (zxlink)
#   zfdefault  -> D_z-f, compact formulation, default CPLEX solver     (zflink)
#   zxgsec     -> D_z-x, Branch-and-Cut with GSEC only, no Benders     (zxlink)
#   zftbend    -> D_z-f, Branch-and-Cut with GSEC + tailored Benders   (zflink)
#   zfautobend -> D_z-f, compact formulation, CPLEX automatic Benders  (zflink)
# zxlink/zflink are derived from MODEL so they can never both be 1 (redundant
# double bound on z) or both be 0 (z left unbounded) as was previously possible
# when they were set independently.
RUN_MODES = {
    "zxdefault":  {"runcode": 0, "model": "zx"},
    "zfdefault":  {"runcode": 0, "model": "zf"},
    "zxgsec":     {"runcode": 2, "model": "zx"},
    "zftbend":    {"runcode": 2, "model": "zf"},
    "zfautobend": {"runcode": 1, "model": "zf"},
}

MODE = "zftbend"  # default when run interactively (i.e. not via main()/CLI)
MODEL = RUN_MODES[MODE]["model"]
zxlink = 1 if MODEL == "zx" else 0
zflink = 1 if MODEL == "zf" else 0




xname = {}
wname = {}
zname = {}
name = {}

EPS = 1e-6
fcuts = 0
ocuts = 0
gtime = 0
ftime = 0
otime = 0
preptime = 0
calltime = 0
total = 0
calls = 0
btight = 0

# The class LazyConsCallback
# allows to add lazy constraints.

def usage(filename1,runcode,capacitym):
    class SECLazyConsCallback(LazyConstraintCallback):

        def cumtimes(self, path):
            """Cumulative travel time from path[0] to every node on path.

            One O(n) sweep replaces the repeated O(n^2) distance() calls that
            previously ran inside the alpha loop (which made that loop O(n^3)).
            """
            dist = self.dist
            cum = {path[0]: 0}
            running = 0
            for a, b in zip(path, path[1:]):
                running += dist[a][b]
                cum[b] = running
            return cum

        def distance(self, path, target):
            """Cumulative time from path[0] to target along path. Now O(n).

            Semantics are unchanged from the original implementation; only the
            repeated path.index() scans have been removed.
            """
            dist = self.dist
            travel = 0
            for a, b in zip(path, path[1:]):
                if a == target:
                    return travel
                travel += dist[a][b]
            return travel



        def tG(self, n1, n2):
            dist = self.dist
            return dist[n1][n2]

        def LB2(self):
            # Depends only on dist/N/Ns/Nd/Nv, all of which are fixed for the
            # instance -- compute once and reuse across every callback call.
            cached = getattr(self, '_lb2_cache', None)
            if cached is not None:
                return cached
            lbval = []
            for j in self.Nd:
                m = self.N.index(j)
                combo = []
                bestsup = {}
                for i in self.Ns:
                    n = self.N.index(i)
                    combo.append([])
                    for k in self.Nv:
                        l = self.N.index(k)
                        combo[i].append(self.dist[l][n] + self.dist[n][m])
                lbval.append(np.min(combo))
            self._lb2_cache = max(lbval)
            return self._lb2_cache


        def linkdistances(self, mode, p):
            # This computation is identical to MILP.linkdistances(); rather
            # than keep two copies (and potentially compute+cache it twice in
            # one run), delegate to the owning MILP instance so it is computed
            # at most once and shared.
            return self.milp.linkdistances(mode, p)

        def __call__(self):
            c0 = time.time()
            p0 = time.time()
            t0 = time.time()
            x = self.x
            w = self.w
            N = self.N
            Ns = self.Ns
            Nd = self.Nd
            Nv = self.Nv
            z = self.z
            vertices = self.vertices

            v = len(x)
            global LPlowerbound
            global fcuts
            global ocuts
            global preptime
            global calltime
            global gtime
            global ftime
            global otime
            global calls
            global btight
            global tight2
            global zflink

            if self.get_num_nodes() == 0:
                LPlowerbound = self.get_best_objective_value()

            # Get the current x solution
            xvij = []
            for a in range(v):
                xvij.append([])
                for b in vertices:
                    xvij[a].append([])
                    for c in vertices:
                        xvij[a][b].append(self.get_values(x[a][b][c]))

            zob = self.get_values(z)

            wvi = []
            route = []
            for a in range(v):
                wvi.append([])
                route.append([])
                for b in vertices:
                    wvi[a].append(self.get_values(w[a][b]))
                    if wvi[a][b] > 0.99:
                        route[a].append(b)

            feasible = []
            round_candidates = []  # (h, path, subtour-node-set) across ALL vehicles this round
            p1 = time.time()
            preptime += p1 - p0

            #################### Prep End

            for h, path in enumerate(route):
                g0 = time.time()
                feasible.append([])
                Gaux2 = nx.DiGraph()
                cnt2 = 0
                for node in path:
                    Gaux2.add_node(N[node], index=cnt2)
                    cnt2 = cnt2 + 1
                for node1 in path:
                    for node2 in path:
                        l1 = N[node1]
                        l2 = N[node2]
                        if xvij[h][node1][node2] > 0.99:
                            Gaux2.add_edge(l1, l2)

                subGraphs2 = nx.weakly_connected_components(Gaux2)  # implemented with bfs

                subtours2 = []
                num_subtours2 = 0

                for c in sorted(subGraphs2, key=len, reverse=True):
                    subtours2.append(c)
                    num_subtours2 += 1
                g1 = time.time()

                gtime += g1 - g0

                if num_subtours2 > 1:
                    feasible[h] = 1
                    # the component holding vehicle h's own depot is the route;
                    # path[0] is the first visited vertex (a supplier), not the depot
                    vnode = Nv[h]
                    for c in subtours2:
                        if len(c) > 1 and vnode not in c:
                            round_candidates.append((h, path, c))
                else:
                    feasible[h] = 0

            # Add a GSEC for every violated subtour this round, largest first.
            if round_candidates:
                f0 = time.time()
                round_candidates.sort(key=lambda item: len(item[2]), reverse=True)
                for h, path, c in round_candidates:
                    # generate cut
                    vars_cut = []
                    bvars_cut = []
                    cvars_cut = []
                    for u in c:
                        i = N.index(u)
                        for v in c:
                            if u != v:
                                j = N.index(v)
                                vars_cut += [x[h][i][j]]
                        bvars_cut += [w[h][i]]
                    tempvars = bvars_cut[:]
                    if multicut == 0:
                        bvars_cut.pop()
                        coeffs_cut = [1.0] * len(vars_cut)
                        bcoeffs_cut = [-1.0] * len(bvars_cut)
                        Mcoeffs_cut = coeffs_cut + bcoeffs_cut
                        Mvars_cut = vars_cut + bvars_cut
                        self.add(constraint=cplex.SparsePair(ind=Mvars_cut, val=Mcoeffs_cut), sense="L",
                                 rhs=0.0)
                        fcuts += 1
                    else:
                        vamod = []
                        for i in tempvars:
                            vamod = [a for a in bvars_cut if a != i]
                            Mvars_cut = vars_cut + vamod
                            coeffs_cut = [1.0] * len(vars_cut)
                            bcoeffs_cut = [-1.0] * len(vamod)
                            Mcoeffs_cut = coeffs_cut + bcoeffs_cut
                            self.add(constraint=cplex.SparsePair(ind=Mvars_cut, val=Mcoeffs_cut), sense="L",
                                     rhs=0.0)
                            fcuts += 1
                f1 = time.time()
                ftime += f1 - f0
            # Tailored Benders optimality cuts only apply to the D_z-f model:
            # when zxlink is active (D_z-x), z is already bounded by the
            # static C13b constraint and GSECs alone are sufficient (paper,
            # "Branch-and-Cut with GSEC"), so this whole alpha/beta pass is
            # skipped for D_z-x runs instead of computing cuts nobody needs.
            if sum(feasible) == 0 and zflink == 1:
                for h, path in enumerate(route):
                    o0 = time.time()
                    feasible[h] = 0
                    successor = {}
                    for i in path:
                        for j in path:
                            if i!=j:
                                out = xvij[h][i][j]
                                if out > 0.1:
                                    successor[i] = j

                    troute = []
                    troute.append(N.index(Nv[h]))
                    for i in path:
                        troute.append(successor[troute[-1]])
                    troute.pop()

                    al = {}
                    alpha = {}
                    beta = {}
                    trace = []
                    for u in troute:
                        trace.append(N[u])

                    # for m in list(set(troute) - set(veh)):
                    # One O(n) sweep gives cumulative arrival times and O(1)
                    # membership/position lookups for the whole alpha loop.
                    cum = self.cumtimes(troute)
                    pos = {node: idx for idx, node in enumerate(troute)}
                    m = troute[len(troute) - 1]
                    al[m] = cum[m]
                    if zob[0] + EPS < cum[m]:
                        alpha_m = cum[m]
                        ind_m = pos[m]
                        for i in range(len(N)):
                            ind_i = pos.get(i)
                            if ind_i is None:
                                alpha[m, i] = alpha_m - self.tG(i, m)
                            else:
                                if ind_m < ind_i:
                                    alt1 = cum[i]
                                    alt2 = alpha_m - self.tG(i, m)
                                    alpha[m, i] = min(alt1, alt2)

                                elif ind_m >= ind_i:
                                    alpha[m, i] = cum[i]

                        if betatight == 1:
                            #print(troute)
                            dep = troute[0]
                            dists = []
                            for i in self.Ns:
                                sup = self.N.index(i)
                                dists.append(self.tG(dep, sup) + self.tG(sup, m))
                            ell = min(dists)
                            tight = al[m] - ell
                        else:
                            tight = 1000

                        for i in range(len(N)):
                            for j in range(len(N)):
                                val = alpha[m, j] - alpha[m, i] - self.tG(i, j)
                                beta[m, i, j] = max(0, min(val,tight))
                                if beta[m, i, j] + EPS < tight:
                                    btight += 1

                        LB3 = self.linkdistances(mode=2, p=1)
                        lbuse = max(LB3, self.LB2())
                        wvars_cut = [w[h][m]]
                        if tight2 == 1:
                            wcoeff_cut = [al[m]-lbuse]
                        else:
                            wcoeff_cut = [al[m]]
                        xvars_cut = []
                        xcoeff_cut = []
                        xrow = x[h]
                        xvals = xvij[h]
                        for i in range(len(N)):
                            for j in range(len(N)):
                                if i != j:
                                    bij = beta[m, i, j]
                                    if SPARSE_CUTS == 1 and bij == 0:
                                        continue
                                    xcoeff_cut += [-bij]
                                    xvars_cut += [xrow[i][j]]
                                    if DEBUG_CUTS == 1 and round(xvals[i][j], 1) != 0 and round(-bij, 1) != 0:
                                        print('i = ' + str(i) + ',' + 'j = ' + str(j) + ':' + 'beta = ' + str(
                                            round(bij, 2)) + ',' + 'xij = ' + str(xvals[i][j]))
                                        print('route', troute)
                                        print('alpha j', str(alpha[m, j]))
                                        print('alpha i', str(alpha[m, i]))
                                        print('alpha m', str(cum[m]))
                                        print('tG i m', str(self.tG(j, m)))
                                        print('====')
                        tvars_cut = wvars_cut + xvars_cut + z
                        tcoeffs_cut = wcoeff_cut + xcoeff_cut + [-1]

                        if DEBUG_CUTS == 1:
                            pcoeff = [round(num, 2) for num in tcoeffs_cut]
                            tbvars = [name[items] for items in tvars_cut]
                            #print(tbvars)
                            #print(pcoeff)
                        # print('current z ',zob)
                        if tight2 == 1:
                            self.add(constraint=cplex.SparsePair(ind=tvars_cut, val=tcoeffs_cut), sense="L", rhs=-lbuse)
                        else:
                            self.add(constraint=cplex.SparsePair(ind=tvars_cut, val=tcoeffs_cut), sense="L", rhs=0.0)
                        ocuts += 1
                        #sys.exit(23)

                o1 = time.time()
                otime += o1 - o0
            c1 = time.time()
            calls += 1
            calltime += c1 - c0

    class CutsFromTable(UserCutCallback):

        # Callback constructor. The fields 'cutlhs' and 'cutrhs' are set externally
        # after the callback is registered.
        def __init__(self, env):
            super().__init__(env)

        def __call__(self):
            for lhs, rhs in zip(self.cutlhs, self.cutrhs):
                # Compute activity of left-hand side
                act = sum(c * x for c, x in zip(lhs.val, self.get_values(lhs.ind)))
                if act > rhs + EPS:
                    print('Adding %s [act = %f]' % (str(lhs), act))
                    self.add(cut=lhs, sense="L", rhs=rhs)

    def draw(G):
        pos = nx.spring_layout(G)
        print(G.nodes)
        # A dictionary with nodes as keys and positions as values.
        nx.draw(G, pos=pos, node_size=210, with_labels=True, font_weight='bold')
        #plt.show()

    class MILP:

        def __init__(self, filename, mode):
            self.problem = cplex.Cplex()
            self.filename = filename
            # self.print_data()
            # chrvrp_io.read_dat_file() parses the instance and builds the
            # full N x N distance matrix once.
            self.s, self.v, self.d, self.dist = chrvrp_io.read_dat_file(filename)

            self.Ns = [i for i in range(0, self.s)]
            self.Nd = [i for i in range(100, self.d + 100)]
            self.Nv = [i for i in range(10, self.v + 10)]
            self.N = self.Ns + self.Nd + self.Nv
            self.Ndot = self.Ns + self.Nd
            self.vehicle_cap = capacitym

            self.x = []
            self.w = []
            self.f = []
            self.z = []

            self.createILP(mode)

            # solve problem
            #self.problem.write('milp2.lp')
            self.solve_and_report(mode)

        #################### File related operations #################
        # get_words/read_dat_file live in chrvrp_io.py; plotter() below calls
        # chrvrp_io.get_words directly.

        def print_data(self):
            print("Debug instance: ")
            for node_i in range(self.nb_nodes):
                for node_j in range(self.nb_nodes):
                    print("Node (%1s, %1s):  weight = %17.5g" %
                          (node_i, node_j, self.arc_cost[node_i][node_j]))
                    # print(node_i, node_j, self.arc_cost[node_i][node_j])

        def linkdistances(self, mode, p):
            # candidate/candidate2 depend only on the distance matrix, so the
            # whole O(|Nd|^2 |Ns| |Nv|) sweep is hoisted out of the callback and
            # cached on first use.
            cached = getattr(self, '_linkdist_cache', None)
            if cached is not None:
                candidate, candidate2 = cached
                if mode == 3:
                    return candidate2[p]
                elif mode == 2:
                    return max(candidate.values())
                return None
            shortp = {}
            secleg = {}
            for k in self.Nd:
                z = self.N.index(k)
                for i in self.Ns:
                    x = self.N.index(i)
                    for j in self.Nd:
                        y = self.N.index(j)
                        if x != z and y != z:
                            if (i, j) in secleg:
                                secleg[i, j] = min(secleg[i, j], self.dist[x][z] + self.dist[z][y])
                            else:
                                secleg[i, j] = self.dist[x][z] + self.dist[z][y]
                        for h in self.Nv:
                            w = self.N.index(h)
                            if x != z and y != z:
                                if (h, j) in shortp:
                                    shortp[h, j] = min(shortp[h, j],
                                                       self.dist[w][x] + self.dist[x][z] + self.dist[z][y])
                                else:
                                    shortp[h, j] = self.dist[w][x] + self.dist[x][z] + self.dist[z][y]

            candidate = {}
            for i in self.Nv:
                candidate[i] = min(shortp[i, j] for j in self.Nd)

            candidate2 = {}
            for i in self.Ns:
                candidate2[i] = min(secleg[i, j] for j in self.Nd)

            self._linkdist_cache = (candidate, candidate2)

            if mode == 3:
                return candidate2[p]
            elif mode == 2:
                return max(candidate.values())

        def addDesignVars(self):
            global name
            global xname
            global wname
            global zname
            # create variables
            for v in range(len(self.Nv)):
                self.x.append([])
                for i in range(len(self.N)):
                    self.x[v].append([])
                    for j in range(len(self.N)):
                        varName2 = "x_" + str(self.Nv[v]) + "_" + str(self.N[i]) + "_" + str(self.N[j])
                        self.x[v][i].append(self.problem.variables.get_num())
                        name[self.problem.variables.get_num()] = varName2
                        self.problem.variables.add(obj=[0.0],
                                            lb=[0.0],
                                            ub=[1.0],
                                            types= ['B'],
                                            names=[varName2])

            for v in range(len(self.Nv)):
                self.w.append([])
                for i in range(len(self.N)):
                    varName3 = "w_" + str(self.Nv[v]) + "_" + str(self.N[i])
                    self.w[v].append(self.problem.variables.get_num())
                    name[self.problem.variables.get_num()] = varName3
                    self.problem.variables.add(obj=[0.0],
                                        lb=[0.0],
                                        ub=[1.0],
                                        types=['B'],
                                        names=[varName3])

            self.z.append(self.problem.variables.get_num())
            name[self.problem.variables.get_num()] = 'z'
            self.problem.variables.add(obj=[1.0],
                              lb=[0.0],
                              ub=[10000000.0],
                              names=['z'])

            # Add the constraints
            # C1 = Incoming arcs = Outgoing arcs
            # sum_[i in N] {x[v][i][j]} - sum_[k in N] {x[v][j][k]} = 0 forall v in V,j in N
            for v in range(len(self.Nv)):
                for j in range(len(self.N)):
                    thelvars = []
                    thelcoeffs = []
                    thervars = []
                    thercoeffs = []
                    for i in range(len(self.N)):
                        if i!=j:
                            thelvars.append(self.x[v][i][j])
                            thelcoeffs.append(1)
                    for k in range(len(self.N)):
                        if k!=j:
                            thervars.append(self.x[v][j][k])
                            thercoeffs.append(-1)
                    thevars = thelvars + thervars
                    thecoeffs = thelcoeffs + thercoeffs
                    constrname = "C2_" + str(self.Nv[v]) + '.' + str(self.N[j])
                    self.problem.linear_constraints.add(names=[constrname],
                                               lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                               senses=["E"],
                                               rhs=[0.0])

            # C2 = Sum of all incoming arcs = visit variable
            # sum_[i in N]{x[v][i][j]} = w[v][j] forall v in V,j in N
            for v in range(len(self.Nv)):
                for j in range(len(self.N)):
                    thevars = []
                    thecoeffs = []
                    for i in range(len(self.N)):
                        if i!=j:
                            thevars.append(self.x[v][i][j])
                            thecoeffs.append(1)
                    thevars = thevars + [self.w[v][j]]
                    thecoeffs = thecoeffs + [-1]
                    constrname = "C2_" + str(self.Nv[v]) + '.' + str(self.N[j])
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["E"],
                                                        rhs=[0.0])

            # C3 = For each customer, there is only one incoming arc
            # sum_[v in Nv] sum_[i in N] {x[v][i][j]} = 1 for all j in Nd
            for j in range(len(self.Nd)):
                thevars = []
                thecoeffs = []
                for v in range(len(self.Nv)):
                    for i in range(len(self.N)):
                        if self.N[i]!=self.Nd[j]:
                            thevars.append(self.x[v][i][j+self.s])
                            thecoeffs.append(1)
                constrname = "C3_" + str(self.Nd[j])
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[1.0])
            """for j in range(len(self.Nd)):
                thevars = []
                thecoeffs = []
                for v in range(len(self.Nv)):
                    thevars.append(self.w[v][j+self.s])
                    thecoeffs.append(1)
                constrname = "C3a_" + str(self.Nd[j])
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[1.0])"""

            #C4 = Each vehicle visits at most one supply point (paper da5);
            # "= 1" with CHRVRP_FORCED_VEHICLES=1
            # sum_[j in Ns] {w[v][j]} <= 1 forall v in V
            for v in (self.Nv):
                thevars = []
                thecoeffs = []
                for j in (self.Ns):
                    q = self.Nv.index(v)
                    thevars.append(self.w[v-10][j])
                    thecoeffs.append(1)
                constrname = "C4_" + str(v)
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E" if FORCED_VEHICLES else "L"],
                                           rhs=[1.0])

            # C5 = Vehicle cannot visit customer directly (needed once C4 is "<= 1")
            # sum_[j in Nd] {x[v][v][j]} = 0 forall v in V
            for v in (self.Nv):
                thevars = []
                thecoeffs = []
                for j in (self.Nd):
                    if v != j:
                        q = self.Nv.index(v)
                        r = self.Nd.index(j)
                        thevars.append(self.x[v-10][len(self.Ndot) + q][len(self.Ns) + r])
                        thecoeffs.append(1)
                constrname = "C5_" + str(v)
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[0.0])

            # C6 = Vehicle cannot go from a customer back to a supply point
            # sum_[i in Ndot] {x[v][i][j]} = 0 forall v in V, j in Ns
            for v in range(len(self.Nv)):
                for j in range(len(self.Ns)):
                    thevars = []
                    thecoeffs = []
                    for i in range(len(self.Ndot)):
                        if i!=j:
                            thevars.append(self.x[v][i][j])
                            thecoeffs.append(1)
                    constrname = "C6_" + str(self.Nv[v]) + '.' + str(self.N[j])
                    self.problem.linear_constraints.add(names=[constrname],
                                               lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                               senses=["E"],
                                               rhs=[0.0])

            # C7 = A vehicle starting from a depot cannot go back to any other depot
            # sum_[i in N] {x[v][i][vdash]} = 0 forall v, vdash in V, v!=vdash
            for v in range(len(self.Nv)):
                for vdash in range(len(self.Nv)):
                    thevars = []
                    thecoeffs = []
                    if v != vdash:
                        for i in range(len(self.N)):
                            q = v + len(self.Ndot)
                            r = vdash + len(self.Ndot)
                            thevars.append(self.x[v][i][r])
                            thecoeffs.append(1)
                        constrname = "C7_" + str(self.Nv[v]) + '.' + str(self.Nv[vdash])
                        self.problem.linear_constraints.add(names=[constrname],
                                                   lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                   senses=["E"],
                                                   rhs=[0.0])
            # C8 = If an arc or a reverse arc exists, then it needs to be visited by that vehicle
            # x[v][i][j] + x[v][j][i] = w[v][i] forall v in V, i,j in N
            for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    for j in range(len(self.N)):
                        thevars = []
                        thecoeffs = []
                        if i!=j:
                            thevars.append(self.x[v][i][j])
                            thecoeffs.append(1)
                            thevars.append(self.x[v][j][i])
                            thecoeffs.append(1)
                            thevars.append(self.w[v][i])
                            thecoeffs.append(-1)
                            constrname = "C8_" + str(self.Nv[v]) + '.' + str(self.N[i])+ '.' + str(self.N[j])
                            self.problem.linear_constraints.add(names=[constrname],
                                                       lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                       senses=["L"],
                                                       rhs=[0.0])

            # C9 = Capacity constraint
            # sum_[i in N] {w[v][i]} + Y+2 forall v in V
            for v in range(len(self.Nv)):
                thevars = []
                thecoeffs = []
                for i in range(len(self.N)):
                    thevars.append(self.w[v][i])
                    thecoeffs.append(1)
                constrname = "C9_" + str(self.Nv[v]) + '.' + str(self.N[i])
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["L"],
                                           rhs=[self.vehicle_cap + 2])
            # C10 = If any node is visted by a vehicle, then the vehicle should have started
            # w[v][i] <= w[v][d(v)] forall v in V, i in N
            for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    thelvars = []
                    thelcoeffs = []
                    thervars = []
                    thercoeffs = []
                    if v+self.s+self.d!= i:
                        thelvars.append(self.w[v][i])
                        thelcoeffs.append(1)
                        thervars.append(self.w[v][v+self.s+self.d])
                        thercoeffs.append(-1)
                        thevars = thelvars + thervars
                        thecoeffs = thelcoeffs + thercoeffs
                        constrname = "C10_" + str(self.Nv[v]) + '.' + str(self.N[i])
                        self.problem.linear_constraints.add(names=[constrname],
                                                   lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                   senses=["L"],
                                                   rhs=[0.0])

            # C13b = Objective function value is the length of the longest route
            # z >= sum_[i,j in N] {x[v][i][j]} forall v in V

            if zxlink == 1:
                for v in range(len(self.Nv)):
                    thevars = []
                    thecoeffs = []
                    for i in range(len(self.N)):
                        for j in self.Ndot:
                            if self.N[i] != j:
                                k = self.N.index(j)
                                thevars.append(self.x[v][i][k])
                                thecoeffs.append(self.dist[i][k])
                    thevars = thevars + self.z
                    thecoeffs = thecoeffs + [-1]
                    constrname = "C13b_" + str(self.Nv[v])
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["L"],
                                                        rhs=[0])


            if lowerb1 == 1:
                print('Lower Bound 1: Active')
                # LB1 (paper): z >= sum_s [w_s t(k_d,s) + (w_m + w_s - 1) t(s,m)], i.e.
                # sum_s w_s (t(k_d,s) + t(s,m)) + w_m sum_s t(s,m) - z <= sum_s t(s,m)
                for v in range(len(self.Nv)):
                    dep = self.s + self.d + v
                    for j in self.Nd:
                        m = self.N.index(j)
                        tsm = [self.dist[i][m] for i in range(len(self.Ns))]
                        thevars = [self.w[v][i] for i in range(len(self.Ns))] + [self.w[v][m]] + self.z
                        thecoeffs = [self.dist[dep][i] + tsm[i] for i in range(len(self.Ns))] + [sum(tsm), -1]
                        constrname = "C16a_" + str(self.Nv[v]) + "." + str(self.N[m])
                        self.problem.linear_constraints.add(names=[constrname],
                                                            lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                            senses=["L"],
                                                            rhs=[sum(tsm)])

            if lowerb2 == 1:
                print('Lower Bound 2: Active')
                lbval = []
                thevars = []
                thecoeffs = []
                for j in self.Nd:
                    m = self.N.index(j)
                    combo = []
                    bestsup = {}
                    for i in self.Ns:
                        n = self.N.index(i)
                        combo.append([])
                        for k in self.Nv:
                            l = self.N.index(k)
                            combo[i].append(self.dist[l][n] + self.dist[n][m])
                    lbval.append(np.min(combo))
                thevars += self.z
                thecoeffs += [1.0]
                constrname = "C17b_"
                self.problem.linear_constraints.add(names=[constrname], lin_expr=[cplex.SparsePair(list(self.z), [1.0])],
                                                    senses=["G"],
                                                    rhs=[max(lbval)])
            # mu = min customers every vehicle must serve, kappa = min vehicles used
            # (paper defs). linkdistances() gives 2-customer paths, so LB3/LB4 use
            # min(mu, 2) customers: valid for any mu, exact for mu = 2.
            mu = max(self.d - (self.v - 1) * self.vehicle_cap, 0)
            kappa = -(-self.d // self.vehicle_cap)
            lam1 = {}       # lam1[v]: min time for vehicle v to reach 1 customer via a supplier
            for v in range(len(self.Nv)):
                dep = self.s + self.d + v
                lam1[v] = min(self.dist[dep][i] + self.dist[i][self.N.index(j)]
                              for i in range(len(self.Ns)) for j in self.Nd)

            if lowerb3 == 1:
                print('Lower Bound 3 active (mu=%d, kappa=%d)' % (mu, kappa))
                if mu >= 2:
                    lbval2b = self.linkdistances(mode=2, p=1)
                elif mu == 1:
                    lbval2b = max(lam1.values())
                else:
                    lbval2b = 0.0
                lbval2b = max(lbval2b, sorted(lam1.values())[kappa - 1])
                self.problem.linear_constraints.add(names=["C17c_"],
                                                    lin_expr=[cplex.SparsePair(list(self.z), [1.0])],
                                                    senses=["G"],
                                                    rhs=[lbval2b])

            if lowerb4 == 1 and mu >= 1:
                print('Lower Bound 4 Active (mu=%d)' % mu)
                # LB4 (paper): z >= sum_s (t(k_d,s) + t_{s,mu}) w^k_s, one row per vehicle
                for v in range(len(self.Nv)):
                    thevars = []
                    thecoeffs = []
                    for i in range(len(self.Ns)):
                        if mu >= 2:
                            lbval3 = self.linkdistances(mode=3, p=self.N[i])
                        else:
                            lbval3 = min(self.dist[i][self.N.index(j)] for j in self.Nd)
                        thevars.append(self.w[v][i])
                        thecoeffs.append(self.dist[self.s + self.d + v][i] + lbval3)
                    thevars = thevars + self.z
                    thecoeffs = thecoeffs + [-1]
                    constrname = "C16b_" + str(self.Nv[v])
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["L"],
                                                        rhs=[0.0])

        def addFlowVars(self):
            # Create variables f(k,i, j)
            # f = []
            for v in range(len(self.Nv)):
                self.f.append([])
                for m in range(len(self.Ndot)):
                    self.f[v].append([])
                    for i in range(len(self.N)):
                        self.f[v][m].append([])
                        for j in range(len(self.N)):
                            varName1 = "flow." + str(self.Nv[v]) + "." + str(self.N[m]) + "." + str(self.N[i]) + "." + str(self.N[j])
                            self.f[v][m][i].append(self.problem.variables.get_num())
                            self.problem.variables.add(obj=[0.0],
                                              lb=[0.0],
                                              ub=[1.0],
                                              names=[varName1])

            # C11 = Flow constraints in subproblem
            for v in range(len(self.Nv)):
                for m in range(len(self.Ndot)):
                    for i in range(len(self.N)):
                        thelvars = []
                        thervars = []
                        thelcoeffs = []
                        thercoeffs = []
                        rhsvar = []
                        rhscoeff = []
                        for j in range(len(self.N)):
                            if i !=j:
                                thelvars.append(self.f[v][m][j][i])
                                thelcoeffs.append(1)
                                thervars.append(self.f[v][m][i][j])
                                thercoeffs.append(-1)

                        if self.N[i] == self.Nv[v]:
                            rhsvar.append(self.w[v][m])
                            rhscoeff.append(1)
                        else:
                            if self.N[i] == self.Ndot[m]:
                                rhsvar.append(self.w[v][m])
                                rhscoeff.append(-1)
                            else:
                                pass
                        thevars = thelvars + thervars + rhsvar
                        thecoeffs = thelcoeffs + thercoeffs + rhscoeff
                        constrname = "C13_" + str(self.Nv[v]) + "." + str(self.Ndot[m]) + "." + str(self.N[i])
                        self.problem.linear_constraints.add(names=[constrname],
                                                            lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                            senses=["E"],
                                                            rhs=[0.0])

            # C12 = Link between flow and arc variables
            for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    for j in range(len(self.N)):
                        if i != j:
                            for m in range(len(self.Ndot)):
                                thevars = []
                                thecoeffs = []
                                thevars.append(self.f[v][m][i][j])
                                thecoeffs.append(1)
                                thevars.append(self.x[v][i][j])
                                thecoeffs.append(-1)
                                constrname = "C12_" + str(self.Nv[v]) + '.' + str(self.N[i])+ '.' + str(self.N[j])
                                self.problem.linear_constraints.add(names=[constrname],
                                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                           senses=["L"],
                                                           rhs=[0.0])

            # C13a = Objective function value is the length of the longest route
            # z >= sum_[i,j in N] {f[v][m][i][j]} forall v in V, m in Ndash

            if zflink == 1:
                for v in range(len(self.Nv)):
                    for m in range(len(self.Ndot)):
                        thevars = []
                        thecoeffs = []
                        for i in range(len(self.N)):
                            for j in range(len(self.N)):
                                if i!= j:
                                    #k = self.N.index(self.Ndot[j])
                                    thevars.append(self.f[v][m][i][j])
                                    thecoeffs.append(self.dist[i][j])
                        thevars = thevars + self.z
                        thecoeffs =  thecoeffs + [-1]
                        constrname = "C14_" + str(self.Nv[v]) + str(m)
                        self.problem.linear_constraints.add(names=[constrname],
                                                   lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                   senses=["L"],
                                                   rhs=[0])


        def createILP(self, mode):
            self.problem.objective.set_sense(self.problem.objective.sense.minimize)
            self.addDesignVars()
            """ mode (runcode) = 0 : Compact formulation, flow variables added   (RUN_MODES: zxdefault/zfdefault)
                mode (runcode) = 1 : Automatic benders, flow variables added     (RUN_MODES: zfautobend)
                mode (runcode) = 2 : No flow variables added, tailored benders   (RUN_MODES: zxgsec/zftbend)"""

            if mode == 0:
                print('### using compact formulation (%s) ###' % MODE)
                self.addFlowVars()
            elif mode == 1:
                print('### using Automatic Benders (%s) ###' % MODE)
                self.addFlowVars()
            elif mode == 2:
                print('### using %s ###' % MODE)
            return

        def triplength(self,inroute):
            lengths = []
            dist = self.dist
            for path in inroute:
                #print(path)
                travel = 0
                ind = len(path) - 1
                # print(ind)
                # print(path[:path.index(target)-1])
                for i in path:
                    if path.index(i) < ind:
                        travel = travel + dist[self.N.index(i)][self.N.index(path[path.index(i) + 1])]
                lengths.append(travel)
            return lengths

        def plotter(self,route,plotr,filename):
            ret = []
            continuation = False
            loc_Vx = []
            loc_Vy = []
            loc_Sx = []
            loc_Sy = []
            loc_Dx = []
            loc_Dy = []
            with open(filename) as f:
                for line in f:
                    words = chrvrp_io.get_words(line)
                    start = list(words[0])[0]
                    if start in ['V', 'S', 'D']:
                        if start == 'V':
                            loc_Vx.append(float(words[1]) / 1000)
                            loc_Vy.append(float(words[2]) / 1000)
                        if start == 'S':
                            loc_Sx.append(float(words[1]) / 1000)
                            loc_Sy.append(float(words[2]) / 1000)
                        if start == 'D':
                            loc_Dx.append(float(words[1]) / 1000)
                            loc_Dy.append(float(words[2]) / 1000)
                    else:
                        pass
                # prepare matrices
                s = self.s
                v = self.v
                d = self.d

                Ns = self.Ns
                Nd = self.Nd
                Nv = self.Nv
                N = self.N
                Ndot = Ns + Nd
                locx = loc_Sx + loc_Dx + loc_Vx
                locy = loc_Sy + loc_Dy + loc_Vy
            plt.scatter(loc_Dx[0:], loc_Dy[0:], c='c', s=200)
            plt.scatter(loc_Sx[0:], loc_Sy[0:], c='r', s=200, marker='^')
            plt.scatter(loc_Vx[0:], loc_Vy[0:], c='y', s=200, marker='s')
            for i in range(len(Nd)):
                k = self.Nd[i]
                plt.annotate('$%d$' % k, (loc_Dx[i], loc_Dy[i]))
                # demand = str(demand1[i])
                # plt.annotate(demand, (loc_Dx[i] - 1.5, loc_Dy[i] + 3), c='r')
            for i in range(len(Ns)):
                plt.annotate('$%d$' % i, (loc_Sx[i], loc_Sy[i]))
                # supply = str(s1[i])
                # plt.annotate(supply, (loc_Sx[i] - 1.5, loc_Sy[i] + 3), c='g')
            for i in range(len(Nv)):
                k = self.Nv[i]
                plt.annotate('$%d$' % k, (loc_Vx[i], loc_Vy[i]))

            troute = []
            for path in range(len(route)):
                troute.append([])
                for i in route[path]:
                    troute[path].append(self.N[i])

            plt.title('Max Trip duration is ' + format(max(self.triplength(troute)), '.2f'))

            title = ""
            for v in range(self.v):
                rndr = numpy.random
                rndg = numpy.random
                rndb = numpy.random
                rndr.seed(v)
                rndg.seed(v)
                rndb.seed(v)
                r = rndr.random()
                b = rndb.random()
                g = rndg.random()
                color = (r, g, b)
                title += format(max(self.triplength([troute[v]])), '.2f') +","
                for m in range(len(Ndot)):
                    for i in range(len(N)):
                        for j in range(len(N)):
                            if (v,m,i,j) in plotr:
                                x1 = locx[i]
                                y1 = locy[i]
                                x2 = locx[j]
                                y2 = locy[j]
                                plt.plot([x1, x2], [y1, y2], c=color, alpha=0.9)
            plt.suptitle(title, size='small')

            plt.axis('scaled')
            plt.xlim([0, 100])
            plt.ylim([0, 100])
            plt.show()
            # outgraph = 'graphmilp.png'
            # plt.savefig(outgraph, bbox_inches='tight')

        def solve_and_report(self, mode):
            global total
            # solve problem
            self.problem.parameters.randomseed.set(10)
            # Defaults for infeasible/no-solution status codes (103/106/108),
            # where the gap/bound queries below are skipped: outputs() (a
            # closure over these locals) would otherwise raise a NameError
            # instead of writing a result row.
            gap = 0
            bound = -1
            nodes = -1
            try:
                # comment out these lines if you want to use the compact flow-based formulation
                if mode == 2:
                    if zflink == 1:
                        print('### using GSEC + Tailored Benders (D_z-f) ###')
                    else:
                        print('### using GSEC only, no Benders cuts (D_z-x) ###')
                    lazyConnectivity = self.problem.register_callback(SECLazyConsCallback)
                    lazyConnectivity.x = self.x
                    lazyConnectivity.w = self.w
                    lazyConnectivity.N = self.N
                    lazyConnectivity.Ns = self.Ns
                    lazyConnectivity.Nd = self.Nd
                    lazyConnectivity.Nv = self.Nv
                    lazyConnectivity.z = self.z
                    lazyConnectivity.dist = self.dist
                    lazyConnectivity.vertices = range(len(self.N))
                    lazyConnectivity.milp = self  # shares linkdistances() cache, see SECLazyConsCallback.linkdistances

                    ### Table Callback

                    if table == 1:
                        usercb = self.problem.register_callback(CutsFromTable)
                        # Generate all disaggregated constraints and put them into a
                        # table that is scanned by the callback.
                        usercb.cutlhs = []
                        usercb.cutrhs = []
                        for v in range(len(self.Nv)):
                            for j in self.Nd:
                                m = self.N.index(j)
                                val = [0]
                                for i in range(len(self.Ns)):
                                    thevars = []
                                    thecoeffs = []
                                    thevars.append(self.w[v][i])
                                    thecoeffs.append(self.dist[i][m])
                                    thevars.append(self.w[v][m])
                                    thecoeffs.append(self.dist[i][self.s + m] + self.dist[self.s + self.d + v][i])
                                    # print(thevars,thecoeffs)
                                    val[0] += self.dist[i][self.s + m]
                                    thevars = thevars + self.z
                                    thecoeffs = thecoeffs + [-1]
                                    usercb.cutlhs.append(cplex.SparsePair(thevars,thecoeffs))
                                    usercb.cutrhs.append(val[0])
                        print('$$$',usercb.cutlhs)
                        print('$$$',usercb.cutrhs)
                    #self.problem.write("tailoredvin.lp")

                elif mode == 1:
                    self.problem.parameters.benders.strategy.set(self.problem.parameters.benders.strategy.values.full)
                self.problem.parameters.timelimit.set(
                    int(os.environ.get('CHRVRP_TIMELIMIT', 3600)))
                # 0 = CPLEX default (all cores); set to match bp_chrvrp.py --threads
                self.problem.parameters.threads.set(
                    int(os.environ.get('CHRVRP_THREADS', 0)))

                to0 = time.time()
                resfile = filename1 + 'log'
                filelog = open(resfile, "w")
                self.problem.set_results_stream(filelog)
                self.problem.set_warning_stream(filelog)
                self.problem.set_error_stream(filelog)
                self.problem.set_log_stream(filelog)
                if os.environ.get('CHRVRP_LPREL') == '1':
                    # LP relaxation only (compact models, runcode 0): write it and stop
                    self.problem.set_problem_type(self.problem.problem_type.LP)
                    self.problem.solve()
                    with open(filename1 + 'lprel', 'w') as f:
                        f.write('%.6f\n%s\n' % (self.problem.solution.get_objective_value(),
                                                self.problem.solution.get_status()))
                    return
                #self.problem.parameters.mip.limits.nodes.set(1)
                #self.problem.solve()
                #self.problem.parameters.mip.limits.nodes.set(999999999999)
                self.problem.solve()
                nodes = self.problem.solution.progress.get_num_nodes_processed()

                if self.problem.solution.get_status() != 103 and self.problem.solution.get_status() != 106 and self.problem.solution.get_status() != 108:
                    gap = self.problem.solution.MIP.get_mip_relative_gap()
                    bound = self.problem.solution.MIP.get_best_objective()
                to1 = time.time()

                total += to1 - to0

            except:
                print(sys.exc_info()[0])
                raise

            #print("Solution status = " + str(self.problem.solution.get_status()))  # 1 optimal, 2 feasible, 3 not feasible
            # the following line prints the corresponding string
            #print(self.problem.solution.status[self.problem.solution.get_status()])

            print("Objective value = " + str(self.problem.solution.get_objective_value()))
            #print("Number of connectivity cuts separated = " + str(CutCount))
            #print("Total separation time = " + str(sepTime))
            #print("Lower bound at the root node = " + str(LPlowerbound))

            # print optimal solution
            """valuesds = self.problem.solution.get_values()
            out2 = valuesds[self.z[0]]
            route = []
            temp = {}

            path = []
            for v in range(len(self.Nv)):
                path.append([])
                for i in range(len(self.N)):
                    out = valuesds[self.w[v][i]]
                    if out > 0.1:
                        path[v].append(i)
                        pass
                        # print('Visit vehicle', str(self.Nv[v]), ' node ', str(self.N[i]))

            successor = {}
            for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    for j in range(len(self.N)):
                        out = valuesds[self.x[v][i][j]]
                        if out > 0.1:
                            successor[self.Nv[v], self.N[i]] = self.N[j]
            troute = []
            for v in self.Nv:
                troute.append([v])
                troute[self.Nv.index(v)].append(successor[v, v])
                for i in path[self.Nv.index(v)]:
                    if self.N[i] not in self.Nv:
                        troute[self.Nv.index(v)].append(successor[v, troute[self.Nv.index(v)][-1]])
                troute[self.Nv.index(v)].pop()"""

            """plotr = {}
            for v in range(len(self.Nv)):
                for m in range(len(self.Ndot)):
                    for i in range(len(self.N)):
                        for j in range(len(self.N)):
                            out = valuesds[self.f[v][m][i][j]]
                            if out > 0.01:
                                plotr[v, m, i, j] = out
                                #print('f_' + str(self.Nv[v]) + '_' + str(self.N[self.N.index(self.Ndot[m])]) + '_' + str(self.N[i]) + '_' + str(self.N[j]))
                                #print('f,' + str(self.Nv[v]) + ',' + str(self.N[self.N.index(self.Ndot[m])]) + ',' + str(self.N[i]) + ',' + str(self.N[j])+ ',' + str(out)+ ',' + str(self.dist[i][j]))"""


            """for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    out = valuesds[self.w[v][i]]
                    if out > 0.1:
                        pass
                        #print('w,' + str(self.Nv[v]) + ',' + str(self.N[i]),',',str(out))
            route = []
            for path in range(len(troute)):
                route.append([])
                for i in troute[path]:
                    route[path].append(self.N.index(i))"""

            #print(self.triplength(troute))
            #print(troute)
            #print(route)

            #print('feasibility cuts count:      ', fcuts)
            #print('optimality cuts count:       ', ocuts)
            #print('beta tightened count:        ', btight)
            #self.plotter(route,plotr, filename = filename1)
            def prtimes():
                print('    Prepwork', round(preptime, 4))
                print('    Graph constr', round(gtime,4))
                print('    FeasCuts', round(ftime,4))
                print('    OptCuts', round(otime,4))
                print('  Total CallBack',round(calltime,4))
                print('Grand Total', round(total,4))
                print('No of times call backs init', calls)
                print(round(preptime, 4),round(gtime,4),round(ftime,4),round(otime,4),round(calltime,4),round(total,4),fcuts,ocuts)
            #prtimes()

            def outputs():
                #print(self.triplength(troute))
                #tot = sum(self.triplength(troute))
                #ast = max(self.triplength(troute))
                #avg = tot / self.v
                #print(mean(cresctime))
                resfile = filename1 + 'rs'
                file = open(resfile, "w")
                file.write('SUMMARY ######### current_val,sol_status,total_time,last_gap,last_bound,nodes \n')
                #file.write(filename1 + ',' + str(round(tot,2)) + ',' + str(round(avg,2)) + ',' + str(round(last,2)) + ',' + str(round(mean(cresctime),2)) + '\n')
                #file.write('Total Path :' + str(tot) + '\n')
                #file.write('Avg Vehicle Path :'+ str(avg) + '\n')
                file.write( str(round(self.problem.solution.get_objective_value(),2)) + '\n')
                file.write(str(self.problem.solution.get_status()) + '\n')
                file.write(str(round(total, 2)) + '\n')
                file.write(str(round(gap, 2)) + '\n')
                file.write(str(round(bound, 2)) + '\n')
                file.write(str(nodes) + '\n')
                #print(last)
                #file.write('Max Res Time :'+ str(last) + '\n')
                #file.write('Avg Res Time :'+ str(mean(cresctime)) + '\n')

            outputs()

    MILP(filename1,runcode)

#filename1 = "instance_size_30_inst_0.atsp"
#usage(filename1,capacitym,runcode)

def main():
    global MODE, MODEL, zxlink, zflink
    if len(sys.argv) != 4 or sys.argv[2].strip().lower() not in RUN_MODES:
        print(len(sys.argv))
        print("Check input")
        print("Usage: python expt1_fast.py <instance.tsp> <mode> <capacity>")
        print("  <mode> is one of: " + ", ".join(RUN_MODES))
    else:
        filename = sys.argv[1]
        MODE = sys.argv[2].strip().lower()
        capacity = int(sys.argv[3])
        # The mode name alone determines both the solve method (runcode) and
        # the model (zx/zf link) -- no separate model argument needed.
        MODEL = RUN_MODES[MODE]["model"]
        zxlink = 1 if MODEL == "zx" else 0
        zflink = 1 if MODEL == "zf" else 0
        usage(filename1=filename, runcode=RUN_MODES[MODE]["runcode"], capacitym=capacity)

if __name__ == "__main__":
    main()