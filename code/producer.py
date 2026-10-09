# producer.py -- Lin-Zhang (arXiv 2609.03784) extension check (floating point + integer counts)
# Split hypersurfaces X_b = {F_0(u_0,v_0)+...+F_m(u_m,v_m)=0} in P^{2m+1}, F_s binary forms of degree d.
# Goal: for every d (incl. d=3..6 excluded by Lin-Zhang's d>=7), find the binary sectors r whose
# RATIONAL closure (Galois orbit) has finite braid monodromy, and the very-general primitive
# middle Hodge rank R(d,m) = dim of type-(m,m) part of the finite-monodromy multi-sectors.
# Conventions follow Lin-Zhang Sec. 6 (reduced Burau, column convention, Hermitian form (6.6)).
import numpy as np, itertools, json, math, sys
from math import gcd

def burau(d, t):
    n = d - 1
    mats = []
    for i in range(1, n + 1):          # sigma_i, i=1..d-1, basis e_1..e_{d-1}
        B = np.eye(n, dtype=complex)
        I = i - 1
        B[I, I] = -t                   # e_i -> -t e_i
        if i > 1: B[I, I - 1] = t       # e_{i-1} -> e_{i-1} + t e_i
        if i < n: B[I, I + 1] = 1       # e_{i+1} -> e_{i+1} + e_i
        mats.append(B)
    return mats

def herm(d, t):
    n = d - 1
    H = np.zeros((n, n), dtype=complex)
    for j in range(n):
        H[j, j] = 2 + t + 1 / t
        if j + 1 < n:
            H[j, j + 1] = -(1 + t)
            H[j + 1, j] = -(1 + 1 / t)
    return H

def inv_vec(d, t):
    return np.array([sum(t ** k for k in range(j)) for j in range(1, d)], dtype=complex)

def quotient(d, t, mats):
    # quotient E_t = V / C v_t ; basis = images of e_2..e_{d-1} (v_1 = 1)
    v = inv_vec(d, t); n = d - 1
    def qcoords(x):  # x in V -> coords in quotient
        return np.array([x[j] - x[0] * v[j] for j in range(1, n)])
    Q = []
    for B in mats:
        M = np.zeros((n - 1, n - 1), dtype=complex)
        for jj, j in enumerate(range(1, n)):
            M[:, jj] = qcoords(B[:, j])
        Q.append(M)
    return Q

def group_order(gens, cap=300000, scale=1e6):
    key = lambda M: tuple(np.round(M.view(float).ravel() * scale).astype(np.int64))
    I = np.eye(gens[0].shape[0], dtype=complex)
    gens = gens + [np.linalg.inv(g) for g in gens]
    seen = {key(I)}; frontier = [I]
    while frontier:
        nf = []
        for A in frontier:
            for g in gens:
                B = g @ A
                k = key(B)
                if k not in seen:
                    seen.add(k); nf.append(B)
                    if len(seen) > cap: return None
        frontier = nf
    return len(seen)

def max_abs_eig(gens, trials=200, length=12, seed=1):
    rng = np.random.default_rng(seed); best = 1.0
    for _ in range(trials):
        A = np.eye(gens[0].shape[0], dtype=complex)
        for _ in range(length):
            g = gens[rng.integers(len(gens))]
            A = A @ (g if rng.random() < .5 else np.linalg.inv(g))
        best = max(best, float(np.max(np.abs(np.linalg.eigvals(A)))))
    return best

def lemma61_failures(d):
    bad = []
    for a in range(1, d):
        n = d // gcd(a, d); q = d // n
        orb = {q * u for u in range(1, n) if gcd(u, n) == 1}
        if orb <= {1, d - 1}: bad.append(a)
    return bad

def galois_orbit(d, r):
    n = d // gcd(r, d); q = d // n
    return sorted({q * u for u in range(1, n) if gcd(u, n) == 1})

out = {"lemma61_failures": {}, "sectors": {}, "R": {}, "hodge_vectors": {}}

# 1. Lemma 6.1 failure set
for d in range(3, 61):
    f = lemma61_failures(d)
    if f: out["lemma61_failures"][d] = f
print("Lemma 6.1 fails (some Galois orbit inside {1,d-1}) for d in", sorted(out["lemma61_failures"]))

# 2./3. per (d,r): invariance, signature, monomial Hodge split, finite group?
finite_sector = {}
for d in range(3, 10):
    z = np.exp(2j * np.pi / d)
    for r in range(1, d):
        t = z ** r
        mats = burau(d, t)
        v = inv_vec(d, t)
        fix_err = max(np.abs(B @ v - v).max() for B in mats)
        rec = {"fix_err": fix_err}
        if abs(t + 1) > 1e-9:
            H = herm(d, t)
            inv_err = max(np.abs(B.conj().T @ H @ B - H).max() for B in mats)
            P = np.eye(d - 1, dtype=complex)[:, 1:]
            Hq = P.conj().T @ H @ P
            ev = np.linalg.eigvalsh((Hq + Hq.conj().T) / 2)
            rec.update(inv_err=inv_err, sig=[int((ev > 1e-9).sum()), int((ev < -1e-9).sum())])
        else:
            rec.update(sig="t=-1 (transvections)")
        # monomial split of W_r: u^a v^b, 0<=a,b<=d-2, e=a+b+2, sector e mod d
        lo = sum(1 for a in range(d - 1) for b in range(d - 1) if a + b + 2 == r)
        hi = sum(1 for a in range(d - 1) for b in range(d - 1) if a + b + 2 == d + r)
        rec["monomial_dims_low_high"] = [lo, hi]
        Q = quotient(d, t, mats)
        if d <= 7:
            rec["group_order_cap300k"] = group_order(Q)
        rec["max_abs_eig"] = max_abs_eig(Q)
        rec["galois_orbit"] = galois_orbit(d, r)
        out["sectors"][f"{d},{r}"] = rec
    # rational finite sectors: r != 0 qualifies iff every Galois conjugate sector is definite and t != -1
    fin = [0]
    for r in range(1, d):
        ok = all(isinstance(out["sectors"][f"{d},{s}"]["sig"], list) and 0 in out["sectors"][f"{d},{s}"]["sig"]
                 for s in galois_orbit(d, r))
        if ok: fin.append(r)
    finite_sector[d] = fin
for d in range(10, 41):
    # beyond d=9 use the closed signature formula (6.10) validated above for d<=9
    fin = [0]
    for r in range(1, d):
        defin = lambda s: (2 * s != d) and (s == 1 or s == d - 1)
        if all(defin(s) for s in galois_orbit(d, r)): fin.append(r)
    finite_sector[d] = fin
print("finite rational sectors:", {d: f for d, f in finite_sector.items() if len(f) > 1})

# 4. very-general Hodge rank by DP over binary factors, plus full primitive Hodge vector
def binary_monomials(d):
    return [(a + b + 2) for a in range(d - 1) for b in range(d - 1)]  # list of e-values

def hodge_counts(d, m, allowed=None):
    # returns dict p -> count of monomials with sum e = d(p+1) (primitive H^{2m-p,p}); if allowed given,
    # restrict every factor's sector (e mod d) to allowed set
    from collections import Counter
    single = Counter(e for e in binary_monomials(d) if allowed is None or (e % d) in allowed)
    dist = Counter({0: 1})
    for _ in range(m + 1):
        nd = Counter()
        for s, c in dist.items():
            for e, c2 in single.items():
                nd[s + e] += c * c2
        dist = nd
    return {p: dist.get(d * (p + 1), 0) for p in range(0, 2 * m + 1)}

def R_formula(d, m):
    if d in (3, 4, 6):
        return sum(math.factorial(m + 1) // (math.factorial(k) ** 2 * math.factorial(m + 1 - 2 * k))
                   * (d - 2) ** (2 * k) * (d - 1) ** (m + 1 - 2 * k) for k in range(0, (m + 1) // 2 + 1))
    return (d - 1) ** (m + 1)

for d in range(3, 13):
    for m in range(1, 5):
        hv = hodge_counts(d, m)
        fin = hodge_counts(d, m, set(finite_sector[d]))
        # all finite-sector monomials must sit in the middle type (pure type check)
        off_middle = sum(c for p, c in fin.items() if p != m)
        Rdp = fin[m]
        out["R"][f"{d},{m}"] = {"R_dp": Rdp, "R_formula": R_formula(d, m), "off_middle_finite": off_middle,
                                "h_mm_prim": hv[m], "b_mid_prim": sum(hv.values())}
        out["hodge_vectors"][f"{d},{m}"] = [hv[p] for p in range(2 * m + 1)]
        if m >= 2 and d <= 9:
            print(f"d={d} m={m}: R_dp={Rdp} R_formula={R_formula(d,m)} off-middle={off_middle} "
                  f"h^(m,m)_prim={hv[m]} b_prim={sum(hv.values())}")

for d, r in [(3, 1), (4, 1), (5, 1), (5, 2), (6, 1), (6, 2), (7, 1)]:
    s = out["sectors"][f"{d},{r}"]
    print(f"d={d} r={r}: sig={s['sig']} lo/hi={s['monomial_dims_low_high']} order={s.get('group_order_cap300k')} "
          f"max|eig|={s['max_abs_eig']:.4f} orbit={s['galois_orbit']}")
json.dump(out, open("producer_out.json", "w"), indent=1, default=str)
