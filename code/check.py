# check.py -- independent checker for producer.py (reads producer_out.json from the current folder)
# (a) exact monodromy group orders over Z[i] / Z[zeta_3] (no floating point), for the sectors the
#     producer flags as finite; (b) Fermat-point check with Shioda's Hodge-character criterion
#     (Galois action on characters; independent of Burau/signature code); (c) Chevalley-Weil
#     eigenspace Hodge numbers of y^d = F(x) vs producer signatures.
import json, math, itertools
from math import gcd
from collections import Counter

# ---------- (a) exact arithmetic ----------
# Z[i]: pair (a,b) = a + b i.  Z[w], w = zeta_3, w^2 = -1 - w: pair (a,b) = a + b w.
def mul_gauss(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def mul_eis(x, y):
    a, b = x; c, d = y   # (a+bw)(c+dw) = ac + (ad+bc) w + bd w^2 = (ac-bd) + (ad+bc-bd) w
    return (a*c - b*d, a*d + b*c - b*d)
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def neg(x): return (-x[0], -x[1])
ZERO, ONE = (0, 0), (1, 0)

def matmul(A, B, mul):
    n = len(A)
    return tuple(tuple(_dot(A[i], [B[k][j] for k in range(n)], mul) for j in range(n)) for i in range(n))
def _dot(row, col, mul):
    s = ZERO
    for a, b in zip(row, col):
        if a != ZERO and b != ZERO: s = add(s, mul(a, b))
    return s

def burau_exact(d, t, mul):
    n = d - 1
    mats = []
    for i in range(1, n + 1):
        B = [[ONE if r == c else ZERO for c in range(n)] for r in range(n)]
        I = i - 1
        B[I][I] = neg(t)
        if i > 1: B[I][I-1] = t
        if i < n: B[I][I+1] = ONE
        mats.append(B)
    # invariant vector v_j = 1 + t + ... + t^{j-1}
    v = []; p = ONE; s = ZERO
    for j in range(n):
        s = add(s, p); v.append(s); p = mul(p, t)
    assert v[0] == ONE
    # verify B v = v exactly
    for B in mats:
        Bv = [_dot(B[r], v, mul) for r in range(n)]
        assert Bv == v, "invariant vector check failed"
    Q = []
    for B in mats:
        M = [[None]*(n-1) for _ in range(n-1)]
        for jj, j in enumerate(range(1, n)):
            col = [B[r][j] for r in range(n)]
            for ii, i2 in enumerate(range(1, n)):
                M[ii][jj] = add(col[i2], neg(mul(col[0], v[i2])))
        Q.append(tuple(tuple(r) for r in M))
    return Q

def closure(gens, mul, cap=400000):
    n = len(gens[0])
    I = tuple(tuple(ONE if r == c else ZERO for c in range(n)) for r in range(n))
    seen = {I}; frontier = [I]
    while frontier:
        nf = []
        for A in frontier:
            for g in gens:   # positive generators suffice for a finite group
                B = matmul(g, A, mul)
                if B not in seen:
                    seen.add(B); nf.append(B)
                    if len(seen) > cap: return None
        frontier = nf
    return len(seen)

res = {}
# d=3: zeta_3 = w -> (0,1); zeta_3^2 = -1-w -> (-1,-1)
# d=4: zeta_4 = i -> (0,1); zeta_4^3 = -i -> (0,-1)   (Gaussian)
# d=6: zeta_6 = 1+w -> (1,1); zeta_6^5 = -w -> (0,-1) (Eisenstein)
cases = [(3, 1, (0, 1), mul_eis), (3, 2, (-1, -1), mul_eis),
         (4, 1, (0, 1), mul_gauss), (4, 3, (0, -1), mul_gauss),
         (6, 1, (1, 1), mul_eis), (6, 5, (0, -1), mul_eis)]
for d, r, t, mul in cases:
    Q = burau_exact(d, t, mul)
    res[f"exact_order d={d} r={r}"] = closure(Q, mul)
    print(f"exact group order d={d} r={r}:", res[f"exact_order d={d} r={r}"])

# ---------- (b) Fermat characters + Shioda criterion ----------
def shioda_hodge(alpha, d):
    N = len(alpha); target = N // 2
    for tt in range(1, d):
        if gcd(tt, d) != 1: continue
        if sum((tt * a) % d for a in alpha) != target * d: return False
    return True

def fermat_check(d, m):
    N = 2 * m + 2
    pairs = [(a, b) for a in range(1, d) for b in range(1, d)]
    bal_total = bal_hodge = lin_total = lin_hodge = 0
    fermat_hodge_total = 0
    for combo in itertools.product(pairs, repeat=m + 1):
        alpha = [x for p in combo for x in p]
        if sum(alpha) % d: continue
        sums = [a + b for a, b in combo]
        hod = shioda_hodge(alpha, d)
        if hod: fermat_hodge_total += 1
        if all(s in (d - 1, d, d + 1) for s in sums) and sums.count(d + 1) == sums.count(d - 1):
            bal_total += 1; bal_hodge += hod
            if all(s == d for s in sums):
                lin_total += 1; lin_hodge += hod
    return dict(balanced_chars=bal_total, balanced_all_hodge=bal_hodge, linear_chars=lin_total,
                linear_hodge=lin_hodge, fermat_hodge_total=fermat_hodge_total)

prod = json.load(open("producer_out.json"))
for d in range(3, 10):
    for m in (2, 3):
        if (d - 1) ** (2 * m + 2) > 3e7: continue
        fc = fermat_check(d, m)
        R = prod["R"][f"{d},{m}"]["R_dp"]
        fc["R_producer"] = R
        res[f"fermat d={d} m={m}"] = fc
        print(f"d={d} m={m}: {fc}")

# ---------- (c) Chevalley-Weil for y^d = prod_{j=1}^d (x - z_j) (d branch points, local exponent 1) ----------
# h^{1,0} of the zeta^{-r} eigenspace = -1 + sum_j <r/d> = r - 1 ; h^{0,1} = d - r - 1
cw_ok = True
for d in range(3, 10):
    for r in range(1, d):
        if 2 * r == d: continue
        sig = prod["sectors"][f"{d},{r}"]["sig"]
        if sorted(sig) != sorted([r - 1, d - r - 1]): cw_ok = False; print("CW mismatch", d, r, sig)
res["chevalley_weil_signature_match_d3to9"] = cw_ok
print("Chevalley-Weil signatures match producer:", cw_ok)
json.dump(res, open("check_out.json", "w"), indent=1)
