"""Aggregated model A_z-f (working copy of
Experiments/Aggregated vs Disagg/NoIndexLB_Working27Aug23.py for the 2026-10 runs).

Changes from the original:
  * experiment switches from the environment, as in expt1_fast.py:
    CHRVRP_LB ("1234" = all lower bounds), CHRVRP_THREADS (0 = all cores),
    CHRVRP_TIMELIMIT (s), CHRVRP_LPREL=1 (solve the LP relaxation only)
  * LB1 implemented as in the paper; LB3/LB4 use mu/kappa (valid when mu < 2)
  * node count added to the .rs file
  * build-time only: instance read once, no milpnoindex.lp dump, flow
    variables/constraints added in bulk (identical model)
Usage: python expt_agg.py <instance.tsp> 0 <capacity>
"""
import os
import cplex
import sys
import networkx as nx
import numpy as np
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from cplex.exceptions import CplexError
from time import perf_counter
import matplotlib.pyplot as plt
import numpy
import time

LPlowerbound = -10 ^ 6
CutCount = 0
sepTime = 0

name = {}
wname = {}
zname = {}
uname = {}

t0 = time.time()
timestamp = {}
tcount = 0
gtime = 0
ftime = 0
otime = 0
preptime = 0
calltime = 0
total = 0
nodes = -1

_lbs =os.environ.get('CHRVRP_LB', '1234')
lowerb1 = int('1' in _lbs)
lowerb2 = int('2' in _lbs)
lowerb3 = int('3' in _lbs)
lowerb4 = int('4' in _lbs)

EPS = 1e-6
fcuts = 0
ocuts = 0
# The class LazyConsCallback
# allows to add lazy constraints.

def usage(filename1,runcode,capacitym):
    class SECLazyConsCallback(LazyConstraintCallback):

        def distance(self, path, target):
            travel = 0
            dist = self.dist
            ind = path.index(target)
            if path.index(target) == 0:
                return 0
            elif path.index(target) == 1:
                return dist[path[0]][path[1]]
            else:
                for i in path:
                    if path.index(i) < ind:
                        travel = travel + dist[i][path[path.index(i)+1]]
                return travel



        def tG(self, n1, n2):
            dist = self.dist
            return dist[n1][n2]

        def __call__(self):
            global t0
            global tcount
            global gtime
            global ftime
            global otime
            global preptime
            global calltime

            tcount += 1

            c0 = time.time()

            xtime = round(time.time() - t0,4)
            timestamp['callstart',tcount] = xtime

            #################### Prep Start
            p0 = time.time()

            x = self.x
            w = self.w
            N = self.N
            Ns = self.Ns
            Nd = self.Nd
            Nv = self.Nv
            z = self.z
            vertices = self.vertices
            v = len(Nv)
            dot = len(N)-v
            global LPlowerbound
            global fcuts
            global ocuts

            start_sep = perf_counter()  # clock()

            if self.get_num_nodes() == 0:
                LPlowerbound = self.get_best_objective_value()

            xij = []
            for a in vertices:
                xij.append([])
                for b in vertices:
                    xij[a].append(self.get_values(x[a][b]))

            zob = self.get_values(z)

            #file.write('Call back triggered, z current = ' + str(numpy.round(zob,2)) +'\n')
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
            p1 = time.time()
            preptime += p1-p0

            sup = []
            demd = []
            veh = []
            for i in range(len(N)):
                if i < len(Ns):
                    sup.append(i)
                elif    i < dot:
                    demd.append(i)
                else:
                    veh.append(i)

            #################### Prep End

            for path in route:
                g0 = time.time()
                feasible.append([])
                Gaux2 = nx.DiGraph()
                cnt2 = 0
                nodecnt2 = len(path)
                for node in path:
                    Gaux2.add_node(node, index=cnt2)
                    cnt2 = cnt2 + 1
                for node1 in path:
                    for node2 in path:
                        l1 = node1
                        l2 = node2
                        if xij[node1][node2] > 0.99:
                            Gaux2.add_edge(l1, l2)

                subGraphs2 = nx.weakly_connected_components(Gaux2)  # implemented with bfs



                subtours2 = []
                num_subtours2 = 0

                for c in sorted(subGraphs2, key=len, reverse=True):
                    subtours2.append(c)
                    num_subtours2 += 1
                g1 = time.time()

                gtime += g1-g0

                print('subtours', subtours2)
                print('nsubtours', num_subtours2)

                if num_subtours2 > 1:
                    f0 = time.time()
                    feasible[route.index(path)] = 1
                    # print("check 1")
                    Scap = subtours2[0]
                    for c in subtours2:
                        #file.write(str(c) + '\n')
                        vnode = path[-1]
                        if vnode not in c:
                            # generate cut
                            for node in c:
                                print(node)
                                dsvar = []
                                dscoeff = []
                                for k in range(len(Ns)):
                                    if k not in Scap:
                                        dsvar.append(x[vnode][k])
                                        dscoeff.append(1)
                                uvar = []
                                ucoeff = []
                                xvar = []
                                xcoeff = []
                                wvar = []
                                wcoeff = []

                                for i in list(set(Scap) & set (demd)):
                                    for j in demd:
                                        if j not in Scap:
                                            xvar.append(x[i][j])
                                            xcoeff.append(1)
                                v1 = vnode-dot
                                wvar.append(w[v1][node])
                                wcoeff.append(-1)
                                Mvars_cut = dsvar + uvar + xvar + wvar
                                Mcoeffs_cut = dscoeff + ucoeff + xcoeff + wcoeff
                                add = []
                                cadd = []
                                for vars in Mvars_cut:
                                    add.append(name[vars])
                                    cadd.append(Mcoeffs_cut[Mvars_cut.index(vars)])
                                print(add)
                                print(cadd)
                                fcuts += 1
                                self.add(constraint=cplex.SparsePair(ind=Mvars_cut, val=Mcoeffs_cut), sense="G", rhs=0.0)

                    f1 = time.time()

                    ftime += f1-f0
                else:
                    veh = []
                    #for path in route:
                    h = route.index(path)
                    feasible[h] = 0
                    """successor = {}
                    for i in path:
                        for j in path:
                            out = xvij[route.index(path)][i][j]
                            if out > 0.1:
                                successor[i] = j
    
                    troute = []
                    troute.append(N.index(Nv[h]))
                    for i in path:
                        troute.append(successor[troute[-1]])
                    troute.pop()
    
    
                    veh= [troute[0]]
    
                    al = {}
                    alpha = {}
                    beta = {}
                    #print(troute,'@@@@@')
                    #for m in list(set(troute) - set(veh)):
                    m = troute[len(troute)-1]
                    al[m] = self.distance(troute, m)
                    last_cust = max(al, key=lambda key: al[key])
    
                    if zob[0]  < self.distance(troute, last_cust) - 0.0000001:
                        for i in range(len(N)):
                            if i not in troute:
                                alpha[m,i] = self.distance(troute,m) - self.tG(i,m)
                            else:
                                ind_i = troute.index(i)
                                ind_m = troute.index(m)
                                if ind_m < ind_i:
                                    alt1 = self.distance(troute,i)
                                    alt2 = self.distance(troute,m) - self.tG(i,m)
                                    alpha[m,i] = min(alt1,alt2)
    
                                elif ind_m >= ind_i:
                                    alpha[m,i] = self.distance(troute,i)
    
                        for i in range(len(N)):
                            for j in range(len(N)):
                                val =  alpha[m,j] - alpha[m,i] - self.tG(i,j)
                                beta[m,i,j] = max(0,val)
    
                        wvars_cut = [w[route.index(path)][last_cust]]
                        wcoeff_cut = [al[last_cust]]
                        xvars_cut = []
                        xcoeff_cut = []
                        for i in range(len(N)):
                            for j in range(len(N)):
                                xcoeff_cut += [-beta[m,i,j]]
                                xvars_cut += [x[route.index(path)][i][j]]
                        tvars_cut = wvars_cut + xvars_cut + z
                        tcoeffs_cut = wcoeff_cut + xcoeff_cut + [-1]
    
                        #print('zob',zob,' route length',self.distance(troute,last_cust))
                        self.add(constraint=cplex.SparsePair(ind=tvars_cut, val=tcoeffs_cut), sense="L", rhs=0.0)
                        #print('8_OPTI Cut',cplex.SparsePair(ind=tvars_cut, val=tcoeffs_cut))
                        ocuts += 1
                    o1 = time.time()
                    otime += o1-o0"""
            c1 = time.time()
            xtime = round(time.time()-t0,4)
            calltime += c1-c0
            timestamp['callend', tcount] = xtime

            if sum(feasible) == 0:
                for path in route:
                    veh = []
                    # for path in route:
                    print('@@@@@@@@',route)

                    for a in vertices:
                        for b in vertices:
                            if xij[a][b] > 0.1:
                                print('Arc ', str(a), ' and ', str(b))

                    for a in range(v):
                        for b in vertices:
                            if wvi[a][b] > 0.1:
                                print('Vehicle ', str(a) , ' visits ', str(b) )

                    h = route.index(path)

                    successor = {}
                    for i in path:
                        for j in path:
                            if i!=j:
                                out = xij[i][j]
                                if out > 0.1:
                                    successor[i] = j

                    print('^^^^^^^^^',successor)
                    troute = []
                    troute.append(N.index(Nv[h]))
                    print('&&&3&&&',successor[troute[len(troute)-1]])
                    for i in path:
                        troute.append(successor[troute[-1]])
                    troute.pop()
                    print('&&&&&&&&&',troute)

                    al = {}
                    alpha = {}
                    beta = {}
                    trace = []
                    for u in troute:
                        trace.append(N[u])
    
                    # for m in list(set(troute) - set(veh)):
                    m = troute[len(troute) - 1]
                    al[m] = self.distance(troute, m)
                    last_cust = max(al, key=lambda key: al[key])

                    if zob[0] + EPS < self.distance(troute, m):
                        print('REached Here -2 $$$$$$$$$')
                        for i in range(len(N)):
                            if i not in troute:
                                alpha[m, i] = self.distance(troute, m) - self.tG(i, m)
                            else:
                                ind_i = troute.index(i)
                                ind_m = troute.index(m)
                                if ind_m < ind_i:
                                    alt1 = self.distance(troute, i)
                                    alt2 = self.distance(troute, m) - self.tG(i, m)
                                    alpha[m, i] = min(alt1, alt2)
    
                                elif ind_m >= ind_i:
                                    alpha[m, i] = self.distance(troute, i)
                        print('REached Here -1 $$$$$$$$$')
                        for i in range(len(N)):
                            for j in range(len(N)):
                                val = alpha[m, j] - alpha[m, i] - self.tG(i, j)
                                beta[m, i, j] = max(0, val)
                        print('REached Here 0 $$$$$$$$$')
    
                        wvars_cut = [w[route.index(path)][last_cust]]
                        wcoeff_cut = [al[last_cust]]
                        xvars_cut = []
                        xcoeff_cut = []
                        print('REached Here 1 $$$$$$$$$')
                        for i in range(len(N)):
                            for j in range(len(N)):
                                if i!= j:
                                    xcoeff_cut += [-beta[m, i, j]]
                                    xvars_cut += [x[i][j]]
                        tvars_cut = wvars_cut + xvars_cut + z
                        print('REached Here 2 $$$$$$$$$')
                        tcoeffs_cut = wcoeff_cut + xcoeff_cut + [-1]
                        #print('current z ',zob)
                        self.add(constraint=cplex.SparsePair(ind=tvars_cut, val=tcoeffs_cut), sense="L", rhs=0.0)
                        ocuts += 1
                ### End Debug Help ###"""

    def draw(G):
        pos = nx.spring_layout(G)
        print(G.nodes)
        # A dictionary with nodes as keys and positions as values.
        nx.draw(G, pos=pos, node_size=210, with_labels=True, font_weight='bold')
        #plt.show()

    class MILP:

        def __init__(self, filename, mode):
            self.problem = cplex.Cplex()
            # self.print_data()
            self.distsv, self.distds, self.distdd, self.dist = self.read_dat_file(filename)[:4]
            self.s = len(self.distsv)
            self.v = len(self.distsv[0])
            self.d = len(self.distdd)

            self.Ns = [i for i in range(0, self.s)]
            self.Nd = [i for i in range(100, self.d + 100)]
            self.Nv = [i for i in range(10, self.v + 10)]
            self.N = self.Ns + self.Nd + self.Nv
            self.Ndot = self.Ns + self.Nd
            self.vehicle_cap = capacitym

            self.x = []
            self.w = []
            self.u = []
            self.f = []
            self.z = []

            self.createILP(mode)

            if os.environ.get('CHRVRP_WRITELP'):    # debug: dump the model
                self.problem.write(os.environ['CHRVRP_WRITELP'])
            # solve problem
            self.solve_and_report(mode)

        #################### File related operations #################
        def get_words(self, line):
            """Return a list of the tokens in line."""
            line = line.replace("\t", " ")
            line = line.replace("\v", " ")
            line = line.replace("\r", " ")
            line = line.replace("\n", " ")
            while line.count("  "):
                line = line.replace("  ", " ")
            line = line.strip()
            return [word + " " for word in line.split(" ")]

        def read_dat_file(self, filename): #closed
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
                    words = self.get_words(line)
                    start = list(words[0])[0]
                    if start in ['V','S','D']:
                        if start == 'V':
                            loc_Vx.append(float(words[1]))
                            loc_Vy.append(float(words[2]))
                        if start == 'S':
                            loc_Sx.append(float(words[1]))
                            loc_Sy.append(float(words[2]))
                        if start == 'D':
                            loc_Dx.append(float(words[1]))
                            loc_Dy.append(float(words[2]))
                    else:
                        pass
                # prepare matrices
                s = len(loc_Sx)
                v = len(loc_Vx)
                d = len(loc_Dx)

                Ns = [i for i in range(0, s)]
                Nd = [i for i in range(100, d + 100)]
                Nv = [i for i in range(10, v + 10)]
                N = Ns + Nd + Nv
                locx = loc_Sx + loc_Dx + loc_Vx
                locy = loc_Sy + loc_Dy + loc_Vy

                distds = []
                for i in Nd:
                    distds.append([])
                    for j in Ns:
                        x1 = float(loc_Dx[i - 100])
                        x2 = float(loc_Sx[j])
                        y1 = float(loc_Dy[i - 100])
                        y2 = float(loc_Sy[j])
                        d1 = (x1 - x2) * (x1 - x2)
                        d2 = (y1 - y2) * (y1 - y2)
                        f = np.sqrt(d1 + d2)/1000
                        distds[i - 100].append(f)
                distsv = []
                for i in Ns:
                    distsv.append([])
                    for j in Nv:
                        x1 = float(loc_Sx[i])
                        x2 = float(loc_Vx[j - 10])
                        y1 = float(loc_Sy[i])
                        y2 = float(loc_Vy[j - 10])
                        d1 = (x1 - x2) * (x1 - x2)
                        d2 = (y1 - y2) * (y1 - y2)
                        f = np.sqrt(d1 + d2)/1000
                        distsv[i].append(f)

                distdd = []
                for i in Nd:
                    distdd.append([])
                    for j in Nd:
                        x1 = float(loc_Dx[i - 100])
                        x2 = float(loc_Dx[j - 100])
                        y1 = float(loc_Dy[i - 100])
                        y2 = float(loc_Dy[j - 100])
                        d1 = (x1 - x2) * (x1 - x2)
                        d2 = (y1 - y2) * (y1 - y2)
                        f = np.sqrt(d1 + d2)/1000
                        distdd[i - 100].append(f)

                dist = []
                for i in range(len(N)):
                    dist.append([])
                    for j in range(len(N)):
                        x1 = float(locx[i])
                        x2 = float(locx[j])
                        y1 = float(locy[i])
                        y2 = float(locy[j])
                        d1 = (x1 - x2) * (x1 - x2)
                        d2 = (y1 - y2) * (y1 - y2)
                        f = np.sqrt(d1 + d2) / 1000
                        dist[i].append(f)
                ret.append(distsv)
                ret.append(distds)
                ret.append(distdd)
                ret.append(dist)
            return ret

        def print_data(self):
            print("Debug instance: ")
            for node_i in range(self.nb_nodes):
                for node_j in range(self.nb_nodes):
                    print("Node (%1s, %1s):  weight = %17.5g" %
                          (node_i, node_j, self.arc_cost[node_i][node_j]))
                    # print(node_i, node_j, self.arc_cost[node_i][node_j])

        def linkdistances(self, mode, p):
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

            if mode == 3:
                return candidate2[p]
            elif mode == 2:
                return max(candidate.values())

        def addDesignVars(self):
            global name
            global wname
            global uname
            global zname

            # create variables


            for i in range(len(self.N)):
                self.x.append([])
                for j in range(len(self.N)):
                    varName2 = "x_" + str(i) + "_" + str(j)
                    self.x[i].append(self.problem.variables.get_num())
                    name[self.problem.variables.get_num()] = varName2
                    if i < self.s and j >= self.s + self.d: #supplier to depot
                        self.problem.variables.add(obj=[0.0],
                                                   types=["B"],
                                                   lb=[0.0],
                                                   ub=[0.0],

                                                   names=[varName2])
                    elif i >= self.s + self.d and j >= self.s and j< self.s + self.d:
                        self.problem.variables.add(obj=[0.0],   #depot directly to customer
                                                   types=["B"],
                                                   lb=[0.0],
                                                   ub=[0.0],

                                                   names=[varName2])
                    elif i <self.s + self.d and j <self.s:
                        self.problem.variables.add(obj=[0.0],   #between suppliers and from customer back to supplier
                                                   types=["B"],
                                                   lb=[0.0],
                                                   ub=[0.0],

                                                   names=[varName2])
                    elif i >= self.s + self.d and j >= self.s + self.d:
                        self.problem.variables.add(obj=[0.0],   #between depots
                                                   types=["B"],
                                                   lb=[0.0],
                                                   ub=[0.0],

                                                   names=[varName2])
                    else:
                        self.problem.variables.add(obj=[0.0],
                                                   types=["B"],
                                          lb=[0.0],
                                          ub=[1.0],

                                          names=[varName2])

            for v in range(len(self.Nv)):
                self.w.append([])
                for i in range(len(self.N)):
                    varName3 = "w_" + str(v) + "_" + str(i)
                    self.w[v].append(self.problem.variables.get_num())
                    name[self.problem.variables.get_num()] = varName3
                    self.problem.variables.add(obj=[0.0],
                                               types=["B"],
                                               lb=[0.0],
                                               ub=[1.0],

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

            for j in range(len(self.N)):
                thelvars = []
                thelcoeffs = []
                thervars = []
                thercoeffs = []
                for i in range(len(self.N)):
                    if i!=j:
                        thelvars.append(self.x[i][j])
                        thelcoeffs.append(1)
                for k in range(len(self.N)):
                    if k!=j:
                        thervars.append(self.x[j][k])
                        thercoeffs.append(-1)
                thevars = thelvars + thervars
                thecoeffs = thelcoeffs + thercoeffs
                constrname = "C1_" + str(self.N[j])
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[0.0])

            # C2 = Sum of all incoming arcs = sum of visit variable
            # sum_[i in N]{x[i][j]} = sum_[v in V]{w[v][j]} j in N
            for j in range(len(self.N)):
                thelvars = []
                thelcoeffs = []
                thervars = []
                thercoeffs = []
                for i in range(len(self.N)):
                    if i!=j:
                        thelvars.append(self.x[i][j])
                        thelcoeffs.append(1)
                for v in range(len(self.Nv)):
                    thervars.append(self.w[v][j])
                    thercoeffs.append(-1)
                thevars = thelvars + thervars
                thecoeffs = thelcoeffs + thercoeffs
                constrname = "C2_" + str(self.N[j])
                self.problem.linear_constraints.add(names=[constrname],
                                                    lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                    senses=["E"],
                                                    rhs=[0.0])

            # C3 = For each customer, there is only one incoming arc
            # sum_[i in N] {x[i][j]} = 1 for all j in Nd
            for j in range(len(self.Nd)):
                thevars = []
                thecoeffs = []
                for i in range(len(self.N)):
                    if self.N[i]!=self.Nd[j]:
                        thevars.append(self.x[i][j+self.s])
                        thecoeffs.append(1)
                constrname = "C3_" + str(self.Nd[j])
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[1.0])

            # C4 = Vehicle can vist max one supplier
            # sum_[j in Nd] {x[v][v][j]} = 0 forall v in V
            for v in range(len(self.Nv)):
                thevars = []
                thecoeffs = []
                for i in range(len(self.Ns)):
                    thevars.append(self.x[v + self.s + self.d][i])
                    thecoeffs.append(1)
                constrname = "C4_" + str(self.Nv[v])
                self.problem.linear_constraints.add(names=[constrname],
                                                    lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                    senses=["L"],
                                                    rhs=[1.0])

            # C5 = Vehicle cannot visit customer directly
            # sum_[j in Nd] {x[v][v][j]} = 0 forall v in V
            active5 = 0
            if active5 ==1:
                for v in range(len(self.Nv)):
                    thevars = []
                    thecoeffs = []
                    for j in range(len(self.Nd)):
                        if self.Nv[v] != self.Nd[j]:
                            thevars.append(self.x[v + self.s + self.d][self.s + j])
                            thecoeffs.append(1)
                    constrname = "C5_" + str(self.Nv[v])
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["E"],
                                                        rhs=[0.0])
            # C6 = Vehicle cannot visit supplier after demand node
            # sum_[j in Nd] {x[v][v][j]} = 0 forall v in V
            active6 = 0
            if active6 ==1:
                for j in range(len(self.Ns)):
                    thevars = []
                    thecoeffs = []
                    for i in range(len(self.Ndot)):
                        if i!=j:
                            thevars.append(self.x[i][j])
                            thecoeffs.append(1)
                    constrname = "C6_" + str(j)
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["E"],
                                                        rhs=[0.0])

            # C7 = If an arc or a reverse arc exists, then it needs to be visited by that vehicle
            # x[v][i][j] + x[v][j][i] = w[v][i] forall v in V, i,j in N
            for i in range(len(self.Nd)):
                for j in range(len(self.Nd)):
                    thelvars = []
                    thelcoeffs = []
                    thervars = []
                    thercoeffs = []
                    if i != j:
                        thelvars.append(self.x[i + self.s][j + self.s])
                        thelcoeffs.append(1)
                        thelvars.append(self.x[j + self.s][i + self.s])
                        thelcoeffs.append(1)
                        for v in range(len(self.Nv)):
                            thervars.append(self.w[v][i + self.s])
                            thercoeffs.append(-1)
                        thevars = thelvars + thervars
                        thecoeffs = thelcoeffs + thercoeffs
                        constrname = "C7_" + str(self.Nv[v]) + '.' + str(self.Nd[i]) + '.' + str(self.Nd[j])
                        self.problem.linear_constraints.add(names=[constrname],
                                                            lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                            senses=["L"],
                                                            rhs=[0.0])

            # C8 = Only one vehicle has to visit a customer once
            # sum_[v in Nv] sum_[i in N] {x[v][i][j]} = 1 for all j in Nd
            for i in range(len(self.Nd)):
                thevars = []
                thecoeffs = []
                for v in range(len(self.Nv)):
                    thevars.append(self.w[v][i + self.s])
                    thecoeffs.append(1)
                constrname = "C8_" + str(self.Nd[i])
                self.problem.linear_constraints.add(names=[constrname],
                                                    lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                    senses=["E"],
                                                    rhs=[1.0])
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

            # C10 = Setting infeasible vehicle arx as zero
            active10 = 0
            if active10 == 1:
                thevars = []
                thecoeffs = []
                for v in range(len(self.Nv)):
                    for vdash in range(len(self.Nv)):
                        if v != vdash:
                            thevars.append(self.x[v+self.s +self.d][vdash+self.s +self.d])
                            thecoeffs.append(1)
                    for i in range(len(self.Nd)):
                        thevars.append(self.x[v + self.s + self.d][i + self.s])
                        thecoeffs.append(1)
                constrname = "C10_"
                self.problem.linear_constraints.add(names=[constrname],
                                           lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                           senses=["E"],
                                           rhs=[0.0])

            # C11 = Setting infeasible supplier arx as zero
            active11 = 0
            if active11 == 1:
                thevars = []
                thecoeffs = []
                for i in range(len(self.Ns)):
                    for idash in range(len(self.Ns)):
                        if i != idash:
                            thevars.append(self.x[i][idash])
                            thecoeffs.append(1)
                    for v in range(len(self.Nv)):
                        thevars.append(self.x[i][v + self.s + self.d])
                        thecoeffs.append(1)
                constrname = "C11_"
                self.problem.linear_constraints.add(names=[constrname],
                                                    lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                    senses=["E"],
                                                    rhs=[0.0])

            #C12 = vehicle cannot visit more than one depot
            active12 = 1
            if active12 == 1:
                for v in range(len(self.Nv)):
                    thevars = []
                    thecoeffs = []
                    for vdash in range(len(self.Nv)):
                        if v!=vdash:
                            thevars.append(self.w[v][vdash+self.s + self.d])
                            thecoeffs.append(1)
                    constrname = "C12_"
                    self.problem.linear_constraints.add(names=[constrname],
                                                        lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                        senses=["E"],
                                                        rhs=[0.0])
            # C12.5 = a vehicle visits a node only if it leaves its own depot (paper ag33)
            # w[v][i] <= w[v][d(v)] forall v in V, i in N
            # This used to fix w[v][d(v)] = 1, which forced every vehicle to drive.
            # Vehicles may now stay idle (w[v][d(v)] = 0), as in the paper (ag29) and
            # in D_z-x/D_z-f; optima are unchanged whenever every vehicle is used.
            active125 = 1
            if active125 == 1:
                for v in range(len(self.Nv)):
                    dep = v + self.s + self.d
                    for i in range(len(self.N)):
                        if i != dep:
                            constrname = "C12b_" + str(self.Nv[v]) + '.' + str(self.N[i])
                            self.problem.linear_constraints.add(names=[constrname],
                                                                lin_expr=[cplex.SparsePair([self.w[v][i], self.w[v][dep]], [1, -1])],
                                                                senses=["L"],
                                                                rhs=[0.0])

            singleton = 1
            if singleton == 1:
                for j in range(len(self.N)):
                    for i in range(len(self.N)):
                        if i!=j:
                            thelvars = []
                            thelcoeffs = []
                            thervars = []
                            thercoeffs = []
                            thelvars.append(self.x[i][j])
                            thelcoeffs.append(1)
                            thervars.append(self.x[j][i])
                            thercoeffs.append(1)
                            thevars = thelvars + thervars
                            thecoeffs = thelcoeffs + thercoeffs
                            constrname = "C16_" + str(self.N[j])
                            self.problem.linear_constraints.add(names=[constrname],
                                                       lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                       senses=["L"],
                                                       rhs=[1.0])
            for v in range(len(self.Nv)):
                for i in range(len(self.N[self.s:])):
                    for j in range(len(self.N[self.s:])):
                        if i!=j:
                            thevars = []
                            thecoeffs = []
                            thevars.append(self.x[i+self.s][j+self.s])
                            thecoeffs.append(1)
                            thevars.append(self.w[v][i+self.s])
                            thecoeffs.append(1)
                            thevars.append(self.w[v][j+self.s])
                            thecoeffs.append(-1)
                            constrname = "C17_" + str(self.Nv[v])+ '_'+str(self.N[i])+'_'+ str(self.N[j])
                            self.problem.linear_constraints.add(names=[constrname],
                                                                lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                                senses=["L"],
                                                                rhs=[1.0])

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
            lam1 = {}
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
                for v in range(len(self.Nv)):
                    for i in range(len(self.Ns)):
                        if mu >= 2:
                            lbval3 = self.linkdistances(mode=3, p=self.N[i])
                        else:
                            lbval3 = min(self.dist[i][self.N.index(j)] for j in self.Nd)
                        thevars = []
                        thecoeffs = []
                        thevars.append(self.w[v][i])
                        thecoeffs.append(self.dist[i][self.s + self.d + v] + lbval3)
                        thevars = thevars + self.z
                        thecoeffs = thecoeffs + [-1]
                        constrname = "C16b_" + str(self.Nv[v]) + str(self.Ns[i])
                        self.problem.linear_constraints.add(names=[constrname],
                                                            lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                            senses=["L"],
                                                            rhs=[0.0])




            """# C13b = Objective function value is the length of the longest route
            # z >= sum_[i,j in N] {x[v][i][j]} forall v in V
            thevars = []
            thecoeffs = []
            for i in range(len(self.N)):
                for j in range(len(self.N)):
                    if i != j:
                        thevars.append(self.x[i][j])
                        thecoeffs.append(self.dist[i][j])
            thevars = thevars + self.z
            thecoeffs = thecoeffs + [-1]
            constrname = "AC13b_" + str(self.Nv[v])
            self.problem.linear_constraints.add(names=[constrname],
                                                lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                senses=["L"],
                                                rhs=[0])"""

        def addFlowVars(self):
            # Create variables f(k,i, j)
            # f = []
            # Variables and rows are built in Python lists and added to CPLEX in
            # one call per block (same variables, rows, order and names as the
            # original one-call-per-row code; only model build time changes).
            nN, nM, nV = len(self.N), len(self.Ndot), len(self.Nv)
            start = self.problem.variables.get_num()
            names = []
            for v in range(nV):
                self.f.append([])
                for m in range(nM):
                    self.f[v].append([])
                    for i in range(nN):
                        base = start + len(names)
                        self.f[v][m].append(list(range(base, base + nN)))
                        names += ["flow." + str(self.Nv[v]) + "." + str(self.N[m]) + "." + str(self.N[i])
                                  + "." + str(self.N[j]) for j in range(nN)]
            self.problem.variables.add(obj=[0.0] * len(names), lb=[0.0] * len(names),
                                       ub=[1.0] * len(names), names=names)

            # Flows never enter a depot (paper ag38b): f[v][m][i][j] = 0 for j in Nv.
            # Without this, two vehicles could swap their depot-to-supplier arcs and
            # route each other's flows through the other depot, since (ag34) does not
            # apply at suppliers.
            depots = range(self.s + self.d, nN)
            self.problem.variables.set_upper_bounds(
                [(self.f[v][m][i][j], 0.0) for v in range(nV) for m in range(nM)
                 for i in range(nN) for j in depots if i != j])

            # C20 = Flow constraints in subproblem
            rows, rnames = [], []
            for v in range(nV):
                for m in range(nM):
                    fvm = self.f[v][m]
                    for i in range(nN):
                        thevars = [fvm[j][i] for j in range(nN) if j != i] + \
                                  [fvm[i][j] for j in range(nN) if j != i]
                        thecoeffs = [1] * (nN - 1) + [-1] * (nN - 1)
                        if self.N[i] == self.Nv[v]:
                            thevars.append(self.w[v][m])
                            thecoeffs.append(1)
                        elif self.N[i] == self.Ndot[m]:
                            thevars.append(self.w[v][m])
                            thecoeffs.append(-1)
                        rows.append(cplex.SparsePair(thevars, thecoeffs))
                        rnames.append("C20_" + str(self.Nv[v]) + "." + str(self.Ndot[m]) + "." + str(self.N[i]))
            self.problem.linear_constraints.add(names=rnames, lin_expr=rows,
                                                senses="E" * len(rows), rhs=[0.0] * len(rows))

            # C11 = Link between flow and arc variables
            rows, rnames = [], []
            for i in range(nN):
                for j in range(nN):
                    if i != j:
                        xij = self.x[i][j]
                        for m in range(nM):
                            name_ijm = "C21_" + str(self.N[i]) + '.' + str(self.N[j]) + '.' + str(self.Ndot[m])
                            for v in range(nV):
                                rows.append(cplex.SparsePair([self.f[v][m][i][j], xij], [1, -1]))
                                rnames.append(name_ijm)
            self.problem.linear_constraints.add(names=rnames, lin_expr=rows,
                                                senses="L" * len(rows), rhs=[0.0] * len(rows))


            # C23 = Link between flow and arc variables
            """for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    for j in range(len(self.N)):
                        if i != j:
                            for m in range(len(self.Ndot)):
                                thevars = []
                                thecoeffs = []
                                thevars.append(self.f[v][m][i][j])
                                thecoeffs.append(1)
                                thevars.append(self.w[v][i])
                                thecoeffs.append(-1)
                                constrname = "C23_" + str(self.Nv[v]) + '.' + str(self.N[i]) + '.' + str(self.N[j])
                                self.problem.linear_constraints.add(names=[constrname],
                                                                    lin_expr=[cplex.SparsePair(thevars, thecoeffs)],
                                                                    senses=["L"],
                                                                    rhs=[0.0])"""

            # C22 = Objective function value is the length of the longest route
            arcs = [(i, j) for i in range(nN) for j in range(nN) if i != j]
            arccost = [self.dist[i][j] for i, j in arcs] + [-1]
            rows, rnames = [], []
            for v in range(nV):
                for m in range(nM):
                    fvm = self.f[v][m]
                    rows.append(cplex.SparsePair([fvm[i][j] for i, j in arcs] + self.z, arccost))
                    rnames.append("C14_" + str(self.Nv[v]) + str(m))
            self.problem.linear_constraints.add(names=rnames, lin_expr=rows,
                                                senses="L" * len(rows), rhs=[0.0] * len(rows))


        def createILP(self, mode):
            self.problem.objective.set_sense(self.problem.objective.sense.minimize)
            self.addDesignVars()
            """ mode = 0 : Compact formulation, flow variables added
                mode = 1 : Automatic benders, flow variables added
                mode = 2 : No flow variables added, tailored benders"""

            if mode == 0:
                #print('### using compact formulation ###')
                self.addFlowVars()
            elif mode == 1:
                #print('### using Automatic Benders ###')
                self.addFlowVars()
            elif mode == 2:
                pass
                #print('### using Tailored Benders ###')
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

        def plotter(self,plotr,filename):
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
                    words = self.get_words(line)
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

            """troute = []
            for path in range(len(route)):
                troute.append([])
                for i in route[path]:
                    troute[path].append(self.N[i])"""

            #plt.title('Max Trip duration is ' + format(max(self.triplength(troute)), '.2f'))

            #title = ""
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
                #title += format(max(self.triplength([troute[v]])), '.2f') + ","
                for m in range(len(Ndot)):
                    for i in range(len(N)):
                        for j in range(len(N)):
                            if (v, m, i, j) in plotr:
                                x1 = locx[i]
                                y1 = locy[i]
                                x2 = locx[j]
                                y2 = locy[j]
                                plt.plot([x1, x2], [y1, y2], c=color, alpha=0.9)
            #plt.suptitle(title, size='small')

            plt.axis('scaled')
            plt.xlim([0, 100])
            plt.ylim([0, 100])
            plt.show()

        def plotter2(self,successor,filename):
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
                    words = self.get_words(line)
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

            """troute = []
            for path in range(len(route)):
                troute.append([])
                for i in route[path]:
                    troute[path].append(self.N[i])"""

            #plt.title('Max Trip duration is ' + format(max(self.triplength(troute)), '.2f'))

            #title = ""
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
            #title += format(max(self.triplength([troute[v]])), '.2f') + ","

            for k in successor:
                i = k
                j = successor[k]
                x1 = locx[i]
                y1 = locy[i]
                x2 = locx[j]
                y2 = locy[j]
                plt.plot([x1, x2], [y1, y2], c=color, alpha=0.9)
                #plt.suptitle(title, size='small')

            plt.axis('scaled')
            plt.xlim([0, 100])
            plt.ylim([0, 100])
            plt.show()

        def solve_and_report(self, mode):
            global total
            # solve problem
            self.problem.parameters.randomseed.set(10)
            try:
                # comment out these lines if you want to use the compact flow-based formulation
                if mode == 2:
                    print('### using GSEC + Tailored Benders ###')
                    lazyConnectivity = self.problem.register_callback(SECLazyConsCallback)
                    lazyConnectivity.x = self.x
                    lazyConnectivity.u = self.u
                    lazyConnectivity.w = self.w
                    lazyConnectivity.N = self.N
                    lazyConnectivity.Ns = self.Ns
                    lazyConnectivity.Nd = self.Nd
                    lazyConnectivity.Nv = self.Nv
                    lazyConnectivity.z = self.z
                    lazyConnectivity.dist = self.dist
                    lazyConnectivity.vertices = range(len(self.N))
                    #self.problem.write("milp2.lp")
                elif mode == 1:
                    self.problem.parameters.benders.strategy.set(self.problem.parameters.benders.strategy.values.full)
                self.problem.parameters.timelimit.set(
                    int(os.environ.get('CHRVRP_TIMELIMIT', 3600)))
                # 0 = CPLEX default (all cores)
                self.problem.parameters.threads.set(
                    int(os.environ.get('CHRVRP_THREADS', 0)))
                resfile = 'agg_'+ str(runcode) + '_' + filename1 + 'log'
                filelog = open(resfile, "w")
                self.problem.set_results_stream(filelog)
                self.problem.set_warning_stream(filelog)
                self.problem.set_error_stream(filelog)
                self.problem.set_log_stream(filelog)
                if os.environ.get('CHRVRP_LPREL') == '1':
                    self.problem.set_problem_type(self.problem.problem_type.LP)
                    self.problem.solve()
                    with open(filename1 + 'lprel', 'w') as f:
                        f.write('%.6f\n%s\n' % (self.problem.solution.get_objective_value(),
                                                self.problem.solution.get_status()))
                    sys.exit(0)
                #self.problem.parameters.mip.limits.nodes.set(1)
                #self.problem.solve()
                #self.problem.parameters.mip.limits.nodes.set(999999999999)
                #self.problem.solve()
                to0 = time.time()
                self.problem.solve()
                to1 = time.time()
                global nodes
                nodes = self.problem.solution.progress.get_num_nodes_processed()

                if self.problem.solution.get_status() != 103 and self.problem.solution.get_status() != 106 and self.problem.solution.get_status() != 108:
                    gap = self.problem.solution.MIP.get_mip_relative_gap()
                    bound = self.problem.solution.MIP.get_best_objective()

                total += to1-to0

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
            valuesds = self.problem.solution.get_values()
            """out2 = valuesds[self.z[0]]
            route = []
            temp = {}

            path = []
            for v in range(len(self.Nv)):
                path.append([])
                for i in range(len(self.N)):
                    out = valuesds[self.w[v][i]]
                    if out > 0.1:
                        path[v].append(i)
                        #pass
                        #print('Visit vehicle', str(self.Nv[v]), ' node ', str(self.N[i]))

            successor = {}

            for i in range(len(self.N)):
                for j in range(len(self.N)):
                    out = valuesds[self.x[i][j]]
                    if out > 0.1:
                        successor[i] = j
                        #print('Arc node ', str(self.N[i]), ' node ', str(self.N[j]))
                        #print(out)
            #print(successor)"""
            """troute = []
            for v in self.Nv:
                troute.append([v])
                troute[self.Nv.index(v)].append(successor[v, v])
                for i in path[self.Nv.index(v)]:
                    if self.N[i] not in self.Nv:
                        troute[self.Nv.index(v)].append(successor[v, troute[self.Nv.index(v)][-1]])
                troute[self.Nv.index(v)].pop()
            # print(troute)"""

            """drg = 1
            if drg==1:
                plotr = {}
                for v in range(len(self.Nv)):
                    for m in range(len(self.Ndot)):
                        for i in range(len(self.N)):
                            for j in range(len(self.N)):
                                out = valuesds[self.f[v][m][i][j]]
                                if out > 0.0001:
                                    plotr[v,m,i,j] = out
                                    #print('f_' + str(self.Nv[v]) + '_' + str(self.N[self.N.index(self.Ndot[m])]) + '_' + str(self.N[i]) + '_' + str(self.N[j]))
                                    #if out!=1:
                                        #print(out)
                self.plotter(plotr,filename = filename1)
            elif drg==2:
                pass
                #self.plotter2(successor, filename=filename1)"""

            """for v in range(len(self.Nv)):
                for i in range(len(self.N)):
                    out = valuesds[self.w[v][i]]
                    if out > 0.1:
                        pass
                        #print('w_' + str(self.Nv[v]) + '_' + str(self.N[i]))
            route = []
            for path in range(len(troute)):
                route.append([])
                for i in troute[path]:
                    route[path].append(self.N.index(i))"""

            #print(self.triplength(troute))
            #print(troute)
            #print(route)


            #file.write(str(self.triplength(troute)))

            #print('feasibility cuts count', fcuts)
            #print('optimality cuts count', ocuts)
            #
            #print('',timestamp)
            #print('times', round(preptime, 4),',',round(gtime,4),',',round(ftime,4),',',round(otime,4),',',round(calltime,4),',',round(total,4))

            def outputs():
                #print(self.triplength(troute))
                #tot = sum(self.triplength(troute))
                #ast = max(self.triplength(troute))
                #avg = tot / self.v
                #print(mean(cresctime))
                resfile = 'agg' + str(filename1[:-4]) + '.rs'
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

    MILP(filename1, runcode)


#filename1 = "instance10S2V3Y4PxLxn25.tsp"
#milp = MILP(filename1,1)
#ile.close()

def main():
    """Set default arguments and parse command line."""
    if len(sys.argv) != 4:
        print(len(sys.argv))
        print("Check input")
    else:
        filename = sys.argv[1]
        mode = int(sys.argv[2])
        capacity = int(sys.argv[3])
        usage(filename1 =  filename, runcode= mode, capacitym= capacity)

if __name__ == "__main__":
    main()