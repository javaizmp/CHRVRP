"""Branch-and-price for the CHR-VRP (min-max rescue time), solved to optimality.

Columns are open vehicle routes  k_d -> s -> j_1 -> ... -> j_q  (q <= Q); the
route time c_r is the arrival time at the last customer (no return arc), which
matches the D_z-x objective (C13b) in expt1_fast.py. Since the supplier only
affects how a route starts, the start cost to the first customer j is
h_k(j) = min_s [t(k_d,s) + t(s,j)] and the supplier is recovered as that argmin.

The min-max objective is handled by a threshold (decision) version: for a
target T only routes with c_r <= T are allowed, and the question becomes
whether every customer can be covered by at most one route per vehicle.

  Master (Phase-1 LP at threshold T, per B&B node):
      min  sum_j y_j
      s.t. sum_{r} a_jr lambda_r + y_j = 1      j in N_d      (pi_j)
           sum_{r in R_k(T)} lambda_r   = 1      k in K        (sigma_k)
      lambda, y >= 0;  R_k(T) = routes of vehicle k with c_r <= T
  (every vehicle has an "empty" route, so the convexity rows are always
  satisfiable.)

  Pricing for vehicle k: elementary route with time <= T and <= Q customers
  maximising sum pi_j  (reduced cost = -sum pi_j - sigma_k), solved by a
  numba-compiled labelling algorithm with time/prize/count/visited(+unreachable)
  dominance; a label-limited heuristic pass runs first, the exact pass only
  when it finds nothing.

  Branching on f_jk = sum_{r in R_k} a_jr lambda_r (customer j served by
  vehicle k). Once every f_jk is integral, every positive route of vehicle k
  covers exactly that vehicle's customer set, so any one of them gives an
  integer solution -- assignment branching alone is complete.

  Outer loop: bisection on T between a proven lower bound LB and the best
  solution value UB. A node is pruned as infeasible when the Lagrangian bound
  z_RMP + sum_k min(0, min_rc_k) > 0. A feasible run lowers UB to the found
  makespan (<= T); an infeasible run raises LB to T.

Two semantics for unused vehicles:
  default           : paper constraint (da5) "<= 1"; a vehicle may stay at its
                      depot, and its empty route costs 0.
  --forced-vehicles : every vehicle drives to a supplier (C4 in expt1_fast.py
                      is "= 1"), so an unused vehicle still costs
                      min_s t(k_d, s) -- reproduces zxgsec's optima exactly.

Parallelism (--threads N): the |K| pricing problems of an iteration share the
same duals and are independent, so they run concurrently in a thread pool (the
numba kernel releases the GIL); columns are then added in vehicle order, so the
column pool is identical to a serial run. The pool IP also gets N CPLEX
threads; the master LP stays single-threaded (small, warm-started dual simplex
does not benefit from threads).

Usage (CPLEX Python):
    python bp_chrvrp.py instance20S2V4Y6PxLxn20.tsp 6 [--timelimit 3600] [--threads 6]
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import cplex
import numpy as np
from numba import njit

import chrvrp_io

FEAS_TOL = 1e-6     # Phase-1 objective at or below this = LP feasible
RC_TOL = 1e-7       # reduced cost must be below -RC_TOL to add a column
# A node is declared LP-infeasible only if its bound exceeds FEAS_TOL + nv*RC_TOL:
# pricing ignores columns with reduced cost in (-RC_TOL, 0), which can make the
# computed bound overestimate the true one by at most nv*RC_TOL.
INT_TOL = 1e-6      # integrality tolerance on f_jk
T_TOL = 1e-9        # route time <= T + T_TOL counts as within threshold


class Instance:
    """Customer-indexed view of a CHR-VRP instance.

    Customers are 0..d-1 (position s+j in the chrvrp_io distance matrix),
    vehicles 0..v-1 (depot at position s+d+k), suppliers 0..s-1.
    """

    def __init__(self, path, Q, forced_vehicles=False):
        s, v, d, dist = chrvrp_io.read_dat_file(path)
        self.path = path
        self.ns, self.nv, self.nd, self.Q = s, v, d, Q
        self.forced_vehicles = forced_vehicles
        # customer-to-customer times
        self.t = [[dist[s + i][s + j] for j in range(d)] for i in range(d)]
        # start cost to first customer j for vehicle k, and the supplier used
        self.h, self.hsup = [], []
        self.empty_cost, self.empty_sup = [], []
        for k in range(v):
            dep = s + d + k
            hk, sk = [], []
            for j in range(d):
                best = min(range(s), key=lambda i: dist[dep][i] + dist[i][s + j])
                hk.append(dist[dep][best] + dist[best][s + j])
                sk.append(best)
            self.h.append(hk)
            self.hsup.append(sk)
            near = min(range(s), key=lambda i: dist[dep][i])
            self.empty_sup.append(near if forced_vehicles else None)
            self.empty_cost.append(dist[dep][near] if forced_vehicles else 0.0)
        if -(-d // Q) > v:
            raise ValueError("infeasible instance: %d customers need more than %d vehicles of capacity %d"
                             % (d, v, Q))
        self.t_arr = np.array(self.t, dtype=np.float64).reshape(d, d)
        self.h_arr = np.array(self.h, dtype=np.float64)

    def route_time(self, k, seq):
        if not seq:
            return self.empty_cost[k]
        c = self.h[k][seq[0]]
        for a, b in zip(seq, seq[1:]):
            c += self.t[a][b]
        return c

    def best_route(self, k, custs):
        """Exact minimum-time route of vehicle k over the customer set custs
        (Held-Karp over the set; |custs| <= Q is small)."""
        custs = list(custs)
        n = len(custs)
        if n == 0:
            return (), self.empty_cost[k]
        if n == 1:
            return (custs[0],), self.h[k][custs[0]]
        t, h = self.t, self.h[k]
        full = (1 << n) - 1
        INF = float("inf")
        f = [[INF] * n for _ in range(1 << n)]
        par = [[-1] * n for _ in range(1 << n)]
        for i in range(n):
            f[1 << i][i] = h[custs[i]]
        for m in range(1, 1 << n):
            fm = f[m]
            for i in range(n):
                fi = fm[i]
                if fi == INF:
                    continue
                ci = custs[i]
                for j in range(n):
                    if m >> j & 1:
                        continue
                    m2 = m | (1 << j)
                    val = fi + t[ci][custs[j]]
                    if val < f[m2][j]:
                        f[m2][j] = val
                        par[m2][j] = i
        last = min(range(n), key=lambda i: f[full][i])
        best = f[full][last]
        order, m, i = [], full, last
        while i != -1:
            order.append(custs[i])
            pi = par[m][i]
            m ^= 1 << i
            i = pi
        return tuple(reversed(order)), best

    def supplier_of(self, k, seq):
        return self.hsup[k][seq[0]] if seq else self.empty_sup[k]


class Route:
    __slots__ = ("k", "seq", "cost", "mask")

    def __init__(self, inst, k, seq):
        self.k = k
        self.seq = tuple(seq)
        self.cost = inst.route_time(k, self.seq)
        m = 0
        for j in self.seq:
            m |= 1 << j
        self.mask = m


# ---------------------------------------------------------------------------
# Heuristic upper bound
# ---------------------------------------------------------------------------

def greedy_solution(inst):
    """Insertion heuristic minimising the makespan, then relocate moves off the
    longest route, then Held-Karp re-optimisation of every route."""
    K, Q = range(inst.nv), inst.Q
    seqs = [[] for _ in K]
    costs = [inst.empty_cost[k] for k in K]
    order = sorted(range(inst.nd), key=lambda j: -min(inst.h[k][j] for k in K))
    for j in order:
        best = None
        for k in K:
            if len(seqs[k]) >= Q:
                continue
            for p in range(len(seqs[k]) + 1):
                cand = seqs[k][:p] + [j] + seqs[k][p:]
                c = inst.route_time(k, cand)
                span = max([c] + [costs[o] for o in K if o != k])
                key = (span, c - costs[k])
                if best is None or key < best[0]:
                    best = (key, k, cand, c)
        _, k, cand, c = best
        seqs[k], costs[k] = cand, c

    def reopt(k):
        seq, c = inst.best_route(k, seqs[k])
        seqs[k], costs[k] = list(seq), c

    for k in K:
        reopt(k)
    improved = True
    while improved:
        improved = False
        kmax = max(K, key=lambda k: costs[k])
        span = costs[kmax]
        best = None
        for j in seqs[kmax]:
            rest = [x for x in seqs[kmax] if x != j]
            _, crest = inst.best_route(kmax, rest)
            for k in K:
                if k == kmax or len(seqs[k]) >= Q:
                    continue
                _, cnew = inst.best_route(k, seqs[k] + [j])
                newspan = max(crest, cnew)
                if newspan < span - 1e-9 and (best is None or newspan < best[0]):
                    best = (newspan, j, k)
        if best is not None:
            _, j, k = best
            seqs[kmax].remove(j)
            seqs[k].append(j)
            reopt(kmax)
            reopt(k)
            improved = True
    return [Route(inst, k, seqs[k]) for k in K]


def improve_solution(inst, routes):
    """Re-optimise each route of an integer solution over its customer set."""
    out = []
    for r in routes:
        seq, _ = inst.best_route(r.k, r.seq)
        out.append(Route(inst, r.k, seq))
    return out


def makespan(routes):
    return max(r.cost for r in routes)


def check_solution(inst, routes):
    seen = []
    assert sorted(r.k for r in routes) == list(range(inst.nv)), "one route per vehicle"
    for r in routes:
        assert len(r.seq) <= inst.Q, "capacity"
        assert abs(r.cost - inst.route_time(r.k, r.seq)) < 1e-9
        seen.extend(r.seq)
    assert sorted(seen) == list(range(inst.nd)), "every customer exactly once"


# ---------------------------------------------------------------------------
# Pricing
# ---------------------------------------------------------------------------

def price(inst, k, pi, sigma_k, T, forbid_mask, label_limit=0, max_cols=50):
    """Elementary capacitated orienteering by labelling.

    Returns (columns, min_rc): columns are customer sequences with reduced cost
    < -RC_TOL (best max_cols of them), min_rc is the most negative reduced cost
    found (exact when label_limit == 0, else heuristic). Extensions that cannot
    beat reduced cost 0 are pruned, so min_rc is exact only below 0.
    """
    d = inst.nd
    allowed = np.array([j for j in range(d) if not (forbid_mask >> j) & 1], dtype=np.int64)
    prize = np.asarray(pi, dtype=np.float64)
    thr = -sigma_k + RC_TOL                 # need prize > thr for rc < -RC_TOL
    pos = np.sort(prize[allowed][prize[allowed] > 0])[::-1]
    topadd = np.zeros(inst.Q + 1)           # topadd[r]: best prize from r more customers
    for r in range(1, inst.Q + 1):
        topadd[r] = topadd[r - 1] + (pos[r - 1] if r - 1 < len(pos) else 0.0)
    node, par, pz = _label(inst.t_arr, inst.h_arr[k], prize, allowed, T + T_TOL, inst.Q,
                           thr, topadd, label_limit, (d + 63) // 64)
    good = np.nonzero(pz[1:] > thr)[0] + 1   # skip the root (the empty route)
    if len(good) == 0:
        return [], 0.0
    good = good[np.argsort(-pz[good], kind="stable")][:max_cols]
    cols = []
    for i in good:
        seq = []
        while i > 0:                        # label 0 is the root (depot + supplier)
            seq.append(int(node[i]))
            i = par[i]
        cols.append(tuple(reversed(seq)))
    return cols, float(-pz[good[0]] - sigma_k)


@njit(cache=True, nogil=True)
def _grow1(a, cap):
    b = np.empty(cap, a.dtype)
    b[:a.shape[0]] = a
    return b


@njit(cache=True, nogil=True)
def _grow2(a, cap):
    b = np.empty((cap, a.shape[1]), a.dtype)
    b[:a.shape[0]] = a
    return b


@njit(cache=True, nogil=True)
def _label(t, h, prize, allowed, tmax, Q, thr, topadd, label_limit, W):
    """Labelling kernel. Label 0 is the root (vehicle at its supplier, time 0);
    extending the root to customer j costs h[j], any other extension t[i][j].
    Labels are created level by level (number of customers), so processing
    them in creation order is a valid label-setting order: a label can only be
    dominated by one with no more customers.
    Masks are W words of 64 bits: vm = visited, um = visited | unreachable.
    Returns (node, parent, prize) for every label created."""
    d = t.shape[0]
    one = np.uint64(1)
    full = ~np.uint64(0)
    cap = 1024
    node = np.empty(cap, np.int64)
    par = np.empty(cap, np.int64)
    tm = np.empty(cap, np.float64)
    pz = np.empty(cap, np.float64)
    nn = np.empty(cap, np.int64)
    alive = np.empty(cap, np.bool_)
    nxt = np.empty(cap, np.int64)
    vm = np.zeros((cap, W), np.uint64)
    um = np.zeros((cap, W), np.uint64)
    head = np.full(d, -1, np.int64)
    n_alive = np.zeros(d, np.int64)
    newv = np.empty(W, np.uint64)
    newu = np.empty(W, np.uint64)

    node[0] = -1
    par[0] = -1
    tm[0] = 0.0
    pz[0] = 0.0
    nn[0] = 0
    alive[0] = True
    L = 1
    q = 0
    while q < L:
        if not alive[q]:
            q += 1
            continue
        if q > 0 and (nn[q] >= Q or pz[q] + topadd[Q - nn[q]] <= thr):
            q += 1
            continue                         # no extension can reach rc < 0
        for a in range(allowed.shape[0]):
            l = allowed[a]
            wl = l >> 6
            bl = one << np.uint64(l & 63)
            if um[q, wl] & bl:
                continue
            time = h[l] if q == 0 else tm[q] + t[node[q], l]
            if time > tmax:
                continue
            p = pz[q] + prize[l]
            n = nn[q] + 1
            for w in range(W):
                newv[w] = vm[q, w]
            newv[wl] |= bl
            if n >= Q:
                for w in range(W):
                    newu[w] = full
            else:
                for w in range(W):
                    newu[w] = newv[w]
                for b in range(allowed.shape[0]):
                    m = allowed[b]
                    if time + t[l, m] > tmax:
                        newu[m >> 6] |= one << np.uint64(m & 63)
            # is the new label dominated by an alive label at l?
            dominated = False
            i = head[l]
            while i != -1:
                if alive[i] and tm[i] <= time + 1e-12 and pz[i] >= p - 1e-12 and nn[i] <= n:
                    sub = True
                    for w in range(W):
                        if um[i, w] & ~newu[w]:
                            sub = False
                            break
                    if sub:
                        dominated = True
                        break
                i = nxt[i]
            if dominated:
                continue
            if label_limit > 0 and n_alive[l] >= label_limit:
                continue
            # kill labels at l that the new one dominates
            i = head[l]
            while i != -1:
                if alive[i] and time <= tm[i] + 1e-12 and p >= pz[i] - 1e-12 and n <= nn[i]:
                    sub = True
                    for w in range(W):
                        if newu[w] & ~um[i, w]:
                            sub = False
                            break
                    if sub:
                        alive[i] = False
                        n_alive[l] -= 1
                i = nxt[i]
            if L == cap:
                cap *= 2
                node = _grow1(node, cap)
                par = _grow1(par, cap)
                tm = _grow1(tm, cap)
                pz = _grow1(pz, cap)
                nn = _grow1(nn, cap)
                alive = _grow1(alive, cap)
                nxt = _grow1(nxt, cap)
                vm = _grow2(vm, cap)
                um = _grow2(um, cap)
            node[L] = l
            par[L] = q
            tm[L] = time
            pz[L] = p
            nn[L] = n
            alive[L] = True
            for w in range(W):
                vm[L, w] = newv[w]
                um[L, w] = newu[w]
            nxt[L] = head[l]
            head[l] = L
            n_alive[l] += 1
            L += 1
        q += 1
    return node[:L], par[:L], pz[:L]


# ---------------------------------------------------------------------------
# Restricted master LP
# ---------------------------------------------------------------------------

class Master:
    def __init__(self, inst):
        self.inst = inst
        lp = cplex.Cplex()
        for stream in (lp.set_results_stream, lp.set_log_stream,
                       lp.set_warning_stream, lp.set_error_stream):
            stream(None)
        lp.objective.set_sense(lp.objective.sense.minimize)
        lp.parameters.threads.set(1)
        d, v = inst.nd, inst.nv
        lp.linear_constraints.add(senses="E" * (d + v), rhs=[1.0] * (d + v))
        self.cov_rows = list(range(d))
        self.conv_rows = list(range(d, d + v))
        lp.variables.add(obj=[1.0] * d, lb=[0.0] * d,
                         columns=[cplex.SparsePair(ind=[j], val=[1.0]) for j in range(d)])
        self.nslack = d
        self.lp = lp
        self.routes = []
        self.keys = {}
        self.ub = []        # current upper bound (0/1) of every route column
        for k in range(v):
            self.add_route(Route(inst, k, ()))

    def add_route(self, r):
        key = (r.k, r.seq)
        if key in self.keys:
            return False
        ind = list(r.seq) + [self.conv_rows[r.k]]
        self.lp.variables.add(obj=[0.0], lb=[0.0], ub=[1.0],
                              columns=[cplex.SparsePair(ind=ind, val=[1.0] * len(ind))])
        self.keys[key] = len(self.routes)
        self.routes.append(r)
        self.ub.append(1.0)
        return True

    def set_active(self, T, forbid):
        """Enable exactly the columns with c_r <= T that respect the node's
        branching decisions (forbid[k] = customers vehicle k may not serve)."""
        changes = []
        for i, r in enumerate(self.routes):
            want = 1.0 if (r.cost <= T + T_TOL and not (r.mask & forbid[r.k])) else 0.0
            if want != self.ub[i]:
                self.ub[i] = want
                changes.append((self.nslack + i, want))
        if changes:
            self.lp.variables.set_upper_bounds(changes)

    def solve(self):
        self.lp.solve()
        sol = self.lp.solution
        if sol.get_status() != sol.status.optimal:
            raise RuntimeError("master LP status %s" % sol.get_status_string())
        duals = sol.get_dual_values()
        x = sol.get_values()
        return (sol.get_objective_value(), duals[:self.inst.nd], duals[self.inst.nd:],
                x[self.nslack:])


# ---------------------------------------------------------------------------
# Branch-and-price feasibility search at a fixed threshold T
# ---------------------------------------------------------------------------

class Solver:
    def __init__(self, inst, timelimit, verbose=True, pool_ip_time=10.0, threads=1):
        self.inst = inst
        self.master = Master(inst)
        self.deadline = time.time() + timelimit
        self.verbose = verbose
        self.pool_ip_time = pool_ip_time
        self.threads = threads
        workers = min(threads, inst.nv)
        self.pool = ThreadPoolExecutor(max_workers=workers) if workers > 1 else None
        self.stats = {"nodes": 0, "lp_solves": 0, "lp_time": 0.0, "pricing_time": 0.0,
                      "pool_ip_time": 0.0, "bp_runs": 0, "threads": threads}
        self._pool_ip_size = -1

    def log(self, msg):
        if self.verbose:
            print(msg)
            sys.stdout.flush()

    def timed_out(self):
        return time.time() > self.deadline

    def column_generation(self, T, forbid):
        """Solve the node LP at threshold T. Returns (status, x) with status in
        {'infeasible', 'feasible', 'timeout'}; x are the route values."""
        inst, m = self.inst, self.master
        m.set_active(T, forbid)
        while True:
            if self.timed_out():
                return "timeout", None
            t0 = time.time()
            obj, pi, sigma, x = m.solve()
            self.stats["lp_time"] += time.time() - t0
            self.stats["lp_solves"] += 1
            if obj <= FEAS_TOL:
                return "feasible", x
            t0 = time.time()
            added, exact = 0, False
            for limit in (8, 0):             # heuristic labelling first, then exact
                rc_sum = 0.0
                job = lambda k: price(inst, k, pi, sigma[k], T, forbid[k], label_limit=limit)
                if self.pool is None:
                    results = map(job, range(inst.nv))
                else:
                    results = self.pool.map(job, range(inst.nv))   # yields in vehicle order
                for k, (seqs, min_rc) in enumerate(results):
                    rc_sum += min_rc
                    for seq in seqs:
                        if m.add_route(Route(inst, k, seq)):
                            added += 1
                if added:
                    break
                exact = True
            self.stats["pricing_time"] += time.time() - t0
            infeas_tol = FEAS_TOL + inst.nv * RC_TOL
            if exact:
                # Lagrangian bound with exact pricing: z_LP >= z_RMP + sum_k min_rc_k
                if obj + rc_sum > infeas_tol:
                    return "infeasible", None
            if not added:
                if obj > infeas_tol:
                    return "infeasible", None   # converged with obj > 0
                return "feasible", x

    def pool_ip(self, T):
        """Min-max integer program over the column pool (routes with c_r <= T).
        Returns the best routes found or None."""
        inst, m = self.inst, self.master
        cand = [r for r in m.routes if r.cost <= T + T_TOL]
        if len(cand) == self._pool_ip_size:
            return None
        self._pool_ip_size = len(cand)
        budget = min(self.pool_ip_time, self.deadline - time.time())
        if budget <= 0:
            return None
        t0 = time.time()
        ip = cplex.Cplex()
        for stream in (ip.set_results_stream, ip.set_log_stream,
                       ip.set_warning_stream, ip.set_error_stream):
            stream(None)
        d, v = inst.nd, inst.nv
        ip.objective.set_sense(ip.objective.sense.minimize)
        # rows: cover (d), convexity (v), z >= sum c_r lambda_r per vehicle (v)
        ip.linear_constraints.add(senses="E" * (d + v) + "L" * v, rhs=[1.0] * (d + v) + [0.0] * v)
        cols = []
        for r in cand:
            ind = list(r.seq) + [d + r.k, d + v + r.k]
            val = [1.0] * (len(r.seq) + 1) + [r.cost]
            cols.append(cplex.SparsePair(ind=ind, val=val))
        ip.variables.add(obj=[0.0] * len(cand), types="B" * len(cand), columns=cols)
        ip.variables.add(obj=[1.0], lb=[0.0],
                         columns=[cplex.SparsePair(ind=[d + v + k for k in range(v)], val=[-1.0] * v)])
        ip.parameters.timelimit.set(budget)
        ip.parameters.threads.set(self.threads)
        ip.solve()
        self.stats["pool_ip_time"] += time.time() - t0
        try:
            vals = ip.solution.get_values()
        except cplex.exceptions.CplexError:
            return None
        return [cand[i] for i in range(len(cand)) if vals[i] > 0.5]

    def integer_solution(self, x):
        """If every f_jk is integral, pick one positive route per vehicle."""
        inst, routes = self.inst, self.master.routes
        f = {}
        for i, val in enumerate(x):
            if val > INT_TOL:
                r = routes[i]
                for j in r.seq:
                    f[j, r.k] = f.get((j, r.k), 0.0) + val
        frac = [(abs(val - 0.5), j, k) for (j, k), val in f.items()
                if INT_TOL < val < 1 - INT_TOL]
        if frac:
            _, j, k = min(frac)
            return None, (j, k)
        chosen = {}
        for i, val in enumerate(x):
            if val > INT_TOL:
                r = routes[i]
                if r.k not in chosen or r.cost < chosen[r.k].cost:
                    chosen[r.k] = r
        return [chosen[k] for k in range(inst.nv)], None

    def feasible_at(self, T, root_only=False):
        """Is there an integer solution with makespan <= T?
        Returns ('feasible', routes) | ('infeasible', None) | ('lp_feasible', None)
        (root_only) | ('timeout', None)."""
        inst = self.inst
        self.stats["bp_runs"] += 1
        root = tuple([0] * inst.nv)
        stack = [root]
        nodes = 0
        while stack:
            if self.timed_out():
                return "timeout", None
            forbid = stack.pop()
            nodes += 1
            self.stats["nodes"] += 1
            status, x = self.column_generation(T, forbid)
            if status == "timeout":
                return "timeout", None
            if status == "infeasible":
                continue
            if nodes == 1:
                pool = self.pool_ip(T)
                if pool is not None and makespan(pool) <= T + T_TOL:
                    return "feasible", pool
                if root_only:
                    return "lp_feasible", None
            sol, (branch) = self.integer_solution(x)
            if sol is not None:
                return "feasible", sol
            j, k = branch
            zero = list(forbid)
            zero[k] |= 1 << j                        # vehicle k does not serve j
            one = [fb | (1 << j) if o != k else fb   # only vehicle k serves j
                   for o, fb in enumerate(forbid)]
            stack.append(tuple(zero))
            stack.append(tuple(one))                 # explored first
        return "infeasible", None

    # -----------------------------------------------------------------------

    def solve(self, tol):
        inst = self.inst
        t_start = time.time()
        # trivial lower bounds: every customer needs some vehicle (LB2-style),
        # and (default semantics) every vehicle drives to a supplier
        lb_cust = max(min(inst.h[k][j] for k in range(inst.nv)) for j in range(inst.nd))
        LB = max([lb_cust] + inst.empty_cost) - 1e-9
        best = improve_solution(inst, greedy_solution(inst))
        check_solution(inst, best)
        UB = makespan(best)
        for r in best:
            self.master.add_route(r)
        self.log("initial bounds: LB=%.4f  UB=%.4f (greedy)" % (LB, UB))

        # Phase A: root LP bound by bisection (cheap; no branching)
        lo, hi = LB, UB
        while hi - lo > tol and not self.timed_out():
            T = 0.5 * (lo + hi)
            status, sol = self.feasible_at(T, root_only=True)
            if status == "timeout":
                break
            if status == "infeasible":
                lo = T
            else:
                hi = T
                if status == "feasible":
                    sol = improve_solution(inst, sol)
                    check_solution(inst, sol)
                    if makespan(sol) < UB:
                        best, UB = sol, makespan(sol)
                        hi = min(hi, UB)
            self.log("  root T=%.4f -> %-11s  [LB_root %.4f, UB %.4f]  cols=%d"
                     % (T, status, lo, UB, len(self.master.routes)))
        root_lb = lo
        LB = max(LB, root_lb)
        self.log("root LP bound %.4f  (UB %.4f, gap %.2f%%)  %.1fs"
                 % (root_lb, UB, 100 * (UB - root_lb) / UB, time.time() - t_start))

        # Phase B: exact bisection with full branch-and-price feasibility runs
        status = "optimal"
        while UB - LB > tol:
            if self.timed_out():
                status = "timeout"
                break
            T = 0.5 * (LB + UB)
            n0 = self.stats["nodes"]
            res, sol = self.feasible_at(T)
            if res == "timeout":
                status = "timeout"
                break
            if res == "feasible":
                sol = improve_solution(inst, sol)
                check_solution(inst, sol)
                if makespan(sol) < UB:
                    best, UB = sol, makespan(sol)
            else:
                LB = T
            self.log("  B&P T=%.4f -> %-10s nodes=%-5d [LB %.4f, UB %.4f]  cols=%d  %.1fs"
                     % (T, res, self.stats["nodes"] - n0, LB, UB, len(self.master.routes),
                        time.time() - t_start))
        check_solution(inst, best)
        return {"status": status, "obj": UB, "lb": LB, "root_lb": root_lb,
                "time": time.time() - t_start, "columns": len(self.master.routes),
                "routes": [{"vehicle": r.k, "supplier": inst.supplier_of(r.k, r.seq),
                            "customers": list(r.seq), "time": r.cost} for r in best],
                **self.stats}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("instance")
    ap.add_argument("capacity", type=int)
    ap.add_argument("--timelimit", type=float, default=3600)
    ap.add_argument("--tol", type=float, default=1e-4,
                    help="absolute optimality tolerance on the makespan")
    ap.add_argument("--forced-vehicles", action="store_true",
                    help="expt1_fast.py semantics (C4 = 1): every vehicle drives to a "
                         "supplier, so an unused vehicle costs its depot-supplier time")
    ap.add_argument("--threads", type=int, default=1,
                    help="pricing threads and pool-IP CPLEX threads; 0 = all cores")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    threads = args.threads if args.threads > 0 else os.cpu_count()

    inst = Instance(args.instance, args.capacity, args.forced_vehicles)
    solver = Solver(inst, args.timelimit, verbose=not args.quiet, threads=threads)
    res = solver.solve(args.tol)
    print("%s  status=%s  obj=%.4f  lb=%.4f  root_lb=%.4f  time=%.2fs  nodes=%d  cols=%d"
          % (args.instance, res["status"], res["obj"], res["lb"], res["root_lb"],
             res["time"], res["nodes"], res["columns"]))
    with open(args.instance + ".bp.json", "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
