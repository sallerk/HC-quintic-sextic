"""dims_gf.py -- dimension checks for special members (generating functions; written independently of i1_dims.py).

Jacobian ring of a split form in k pairs: Hilbert series with sector grading
    H(t; z) = prod_s  sum_{e=0}^{2d-4} min(e+1, 2d-3-e) t^e z_s^{(e+2) mod d}.
The coefficient of t^{(p+1)d-2k} z^c is dim V_c^{n-p,p}_prim (n = 2k-2), by Griffiths' theorem and the sector
label c_s = e_s + 2 (mod d) of the note.

(B) identity behind the dimension formula of Theorem B (Theorem 8.4) of the note: for T = {first j indices}, c with c_s = 0 off T:
      h^{p}(V_c(X_F)) = h^{p-(k-j)}(V_{c_T}(X_T)) * (d-1)^{k-j}.
(PSI) Lemma 9.1 of the note: for c with all c_s != 0 (sum 0 mod d):
      h^{p}(V_c) = coefficient of y^{p - (sigma-1)} in prod_s ((c_s-1) + (d-1-c_s) y),
      sigma = sum_s c_s / d (c_s in [1,d-1]); i.e. V_c ~ (x)_s P_s(c_s) (x) (Fermat eigenline of type (k-1-sigma, sigma-1)).
      Also dim V_c = (d-2)^k.
Ranges: 3 <= d <= 9, 1 <= k-1 = m <= 4 (k = m+1 pairs); for (B) also j = 1..k.
"""
import itertools, json, sys
from collections import defaultdict


def series(d, k):
    """dict: (e_total, c tuple) -> dim."""
    one = [(e, min(e + 1, 2 * d - 3 - e), (e + 2) % d) for e in range(0, 2 * d - 3)]
    S = {(0, ()): 1}
    for _ in range(k):
        T = defaultdict(int)
        for (e, c), v in S.items():
            for (e1, h, lab) in one:
                if h > 0:
                    T[(e + e1, c + (lab,))] += v * h
        S = dict(T)
    return S


def hodge_from_series(d, k, S):
    """dict c -> {p: dim} (q-index p of H^{n-p,p})."""
    out = defaultdict(dict)
    for (e, c), v in S.items():
        if (e + 2 * k) % d == 0:
            p = (e + 2 * k) // d - 1
            out[c][p] = out[c].get(p, 0) + v
    return out


def main():
    res = {"B_fail": 0, "B_checked": 0, "PSI_fail": 0, "PSI_checked": 0}
    for d in range(3, 10):
        cache = {}
        for k in range(1, 6):
            cache[k] = hodge_from_series(d, k, series(d, k))
        for k in range(2, 6):  # m = k-1 in 1..4
            HF = cache[k]
            # (B)
            for j in range(1, k + 1):
                HT = cache[j]
                mult = (d - 1) ** (k - j)
                for cT, hT in HT.items():
                    c = cT + (0,) * (k - j)
                    exp = {p + (k - j): v * mult for p, v in hT.items()}
                    got = {p: v for p, v in HF.get(c, {}).items() if v}
                    res["B_checked"] += 1
                    if got != exp:
                        res["B_fail"] += 1
                        if res["B_fail"] < 5:
                            print("B FAIL", d, k, j, cT, got, exp)
            # (PSI)
            for c, h in HF.items():
                if 0 in c:
                    continue
                sigma = sum(c) // d
                poly = {0: 1}
                for cs in c:
                    new = defaultdict(int)
                    for a, v in poly.items():
                        new[a] += v * (cs - 1)
                        new[a + 1] += v * (d - 1 - cs)
                    poly = new
                exp = {a + sigma - 1: v for a, v in poly.items() if v}
                got = {p: v for p, v in h.items() if v}
                res["PSI_checked"] += 1
                if got != exp or sum(got.values()) != (d - 2) ** k:
                    res["PSI_fail"] += 1
                    if res["PSI_fail"] < 5:
                        print("PSI FAIL", d, k, c, got, exp)
    print(json.dumps(res))
    with open(sys.argv[1], "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
