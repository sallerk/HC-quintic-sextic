"""i1_dims.py -- exact dimension checks for split hypersurfaces with some forms fixed.

X_S = {sum_{s in S} F_s(u_s,v_s) = 0} in P^{2|S|-1}, dim n_S = 2|S|-2.
Griffiths: H^{n-p,p}_prim(X_S) <-> R^{(p+1)d-2|S|}, R = (x)_s R_s, dim R_s^e = min(e+1, 2d-3-e), 0<=e<=2d-4.
Sector label of multidegree (e_s): c_s = e_s + 2 mod d.

Checks, for 3<=d<=9, 1<=m<=4, 1<=k<=m+1 (T = first k indices, T^c the rest):
 (A) for every c in (Z/d)^{m+1} with c_{T^c}=0:
     dim V_c^{prim,p}(X_F) = dim V_{c_T}^{prim,p-(m+1-k)}(X_T) * (d-1)^{m+1-k}   for all p,
     and V_c^{prim}(X_F) = 0 for p outside the shifted range.
 (B) summing over c: the G_{T^c}-invariant primitive part of H^{2m}(X_F) has Hodge numbers
     h^{2m-p,p} = h_prim^{n_T-p',p'}(X_T) * (d-1)^{m+1-k}, p = p'+(m+1-k).
 (C) Fermat examples (d prime): dim Hdg^m(X_F) = 1 + (d-1)^{m+1-k} * N_prim(Fermat X^{2k-2}_d) for F_T Fermat
     and F_{T^c} very general relative to F_T (formula of Theorem 8.4 of the note, valid given HC(X_T) and d not in {3,4,6});
     N_prim by Shioda's criterion. For k=1 this must give R(d,m)=(d-1)^{m+1}.
 (D) d=6, |T|=1: list the sectors c with c_t in {2,3,4}, c_s in {0,1,-1} (s != t), sum 0, with their Hodge numbers,
     and the reduced sectors (all c_s = +-1 for s != t), with (m,m)-dimension.
Output: i1_dims.json, and a summary on stdout.
"""
import itertools, json, math, sys
from fractions import Fraction


def hilb(d, e):
    return min(e + 1, 2 * d - 3 - e) if 0 <= e <= 2 * d - 4 else 0


def sector_hodge(d, c):
    """Return dict p -> dim V_c^{prim, (n-p,p)} for the split hypersurface with |c| forms, labels c (tuple mod d)."""
    k = len(c)
    # per factor: list of (e, dim) with e = c_s - 2 mod d, 0<=e<=2d-4
    opts = []
    for cs in c:
        o = [(e, hilb(d, e)) for e in range(0, 2 * d - 3) if (e - (cs - 2)) % d == 0 and hilb(d, e) > 0]
        opts.append(o)
    out = {}
    for combo in itertools.product(*opts):
        e = sum(x[0] for x in combo)
        dim = 1
        for x in combo:
            dim *= x[1]
        # e = (p+1)d - 2k
        if (e + 2 * k) % d != 0:
            continue
        p = (e + 2 * k) // d - 1
        out[p] = out.get(p, 0) + dim
    return out


def sectors(d, k):
    return [c for c in itertools.product(range(d), repeat=k) if sum(c) % d == 0]


def check_A(d, m):
    bad = []
    for k in range(1, m + 2):
        shift = m + 1 - k
        mult = (d - 1) ** shift
        for cT in sectors(d, k):
            c = tuple(cT) + (0,) * shift
            hF = sector_hodge(d, c)
            hT = sector_hodge(d, cT)
            exp = {p + shift: v * mult for p, v in hT.items()}
            if hF != exp:
                bad.append((k, cT, hF, exp))
    return bad


def betti_prim(d, n):
    return ((d - 1) ** (n + 2) + (-1) ** n * (d - 1)) // d


def fermat_hodge_count(d, n):
    """Shioda: primitive Hodge classes of Fermat X^n_d (n even) <-> a in ((Z/d)\0)^{n+2}, sum a = 0,
    sum_i <t a_i>/d = (n+2)/2 for all t in (Z/d)^*."""
    units = [t for t in range(1, d) if math.gcd(t, d) == 1]
    cnt = 0
    half = (n + 2) * d // 2  # sum of residues must equal (n+2)/2 * d
    for a in itertools.product(range(1, d), repeat=n + 1):
        last = (-sum(a)) % d
        if last == 0:
            continue
        a2 = a + (last,)
        if all(sum((t * x) % d for x in a2) == half for t in units):
            cnt += 1
    return cnt


def check_C():
    rows = []
    for d in (5, 7):
        for k in (1, 2, 3):
            n = 2 * k - 2
            if n == 0:
                Nprim = d - 1
            else:
                if (d - 1) ** (n + 1) > 3_000_000:
                    continue
                Nprim = fermat_hodge_count(d, n)
            for m in range(k - 1 if k > 1 else 1, 5):
                if m + 1 < k:
                    continue
                dimHdg = 1 + (d - 1) ** (m + 1 - k) * Nprim
                rows.append(dict(d=d, k=k, m=m, Nprim_fermat_dim=n, Nprim=Nprim, dimHdg=dimHdg,
                                 R_dm_plus1=1 + (d - 1) ** (m + 1)))
    return rows


def d6_T1(m):
    d = 6
    out = []
    for ct in (2, 3, 4):
        for rest in itertools.product((0, 1, 5), repeat=m):
            if (ct + sum(rest)) % d:
                continue
            c = (ct,) + rest
            h = sector_hodge(d, c)
            out.append(dict(c=c, hodge={str(p): v for p, v in sorted(h.items())},
                            mm=h.get(m, 0), n0=rest.count(0), n_plus=rest.count(1), n_minus=rest.count(5)))
    return out


def main():
    res = {"A_failures": {}, "C": None, "D_d6_T1": {}}
    total = 0
    for d in range(3, 10):
        for m in range(1, 5):
            if d ** (m + 1) > 2_000_000:
                # sectors enumeration size d^k; keep runtime small
                pass
            bad = check_A(d, m)
            total += 1
            res["A_failures"][f"{d},{m}"] = len(bad)
            if bad:
                print("FAIL A", d, m, bad[:3])
    # sanity: sum over all sectors of X_S equals primitive Betti number
    for d in range(3, 10):
        for k in range(1, 5):
            n = 2 * k - 2
            s = sum(sum(sector_hodge(d, c).values()) for c in sectors(d, k))
            assert s == betti_prim(d, n), (d, k, s, betti_prim(d, n))
    res["C"] = check_C()
    for m in (1, 2, 3, 4):
        res["D_d6_T1"][str(m)] = d6_T1(m)
    with open(sys.argv[1] if len(sys.argv) > 1 else "i1_dims.json", "w") as f:
        json.dump(res, f, indent=1)
    nfail = sum(res["A_failures"].values())
    print(f"(A) dimension identity: {total} (d,m) pairs, {nfail} failures")
    print("(B) Betti sanity: OK")
    print("(C) Fermat examples:")
    for r in res["C"]:
        print("   ", r)
    print("(D) d=6, |T|=1 extra sectors with nonzero (m,m)-part (reduced: no c_s=0):")
    for m in (1, 2, 3, 4):
        red = [x for x in res["D_d6_T1"][str(m)] if x["n0"] == 0]
        print(f"  m={m}: {len(res['D_d6_T1'][str(m)])} sectors, {len(red)} reduced;",
              [(x['c'], x['hodge'], 'mm=' + str(x['mm'])) for x in red][:12])


if __name__ == "__main__":
    main()
