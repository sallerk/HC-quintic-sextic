"""fermat_T.py -- dim Hdg^m(X_F) when the forms F_T are Fermat (u^d+v^d) and the others are very general relative to them.

The sectors that can carry Hodge classes (Proposition 8.3 of the note) form a local system with finite monodromy over
{F_T fixed} x U_d^{T^c} (P_s(c_s), c_s in Phi_d, s not in T, have finite monodromy), hence (theorem of the fixed
part) a constant Hodge structure along it; so their Hodge classes can be counted at the all-Fermat member.

Method 1 (model): d = 5: 1 + 4^{m+1-j} * N_prim(Fermat X^{2j-2}_5)       (Theorem B, T = j Fermat forms)
                 d = 6, |T| = 1: 1 + R(6,m) + extra, extra from c_t in {2,3,4}: V_c ~ P_t(c_t) (x) (x)_{s!=t} P_s(c_s)
                 (x) L_c, with P_t(g) = sum of eigenlines (a_u,a_v,a_y) of the Fermat sextic curve (a_u+a_v = g),
                 type (1,0) if a_u+a_v+a_y = 6 else (0,1); P(1): 4 x (0,1); P(5): 4 x (1,0); P(0): 5 x (0,0);
                 L_c of type (m - sigma, n0 + sigma - 1), sigma = sum c_s / 6.
Method 2 (Shioda): count characters a in ((Z/d)\0)^{2m+2}, sum 0, with sum_i <u a_i> = (m+1) d for all units u,
                 whose sector c_s = a_{2s} + a_{2s+1} lies in the allowed set; plus 1 for h^m.
"""
import itertools, json, math, sys


def units(d):
    return [u for u in range(1, d) if math.gcd(u, d) == 1]


def shioda_count(d, n, allowed=None):
    """Primitive Hodge classes of Fermat X^n_d (n even), optionally restricted by allowed(c) on pair-sectors."""
    U = units(d)
    half = (n + 2) * d // 2
    cnt = 0
    for a in itertools.product(range(1, d), repeat=n + 1):
        last = (-sum(a)) % d
        if last == 0:
            continue
        a2 = a + (last,)
        if allowed is not None:
            c = tuple((a2[2 * s] + a2[2 * s + 1]) % d for s in range((n + 2) // 2))
            if not allowed(c):
                continue
        if all(sum((u * x) % d for x in a2) == half for u in U):
            cnt += 1
    return cnt


def R6(m):
    return sum(math.factorial(m + 1) // (math.factorial(k) ** 2 * math.factorial(m + 1 - 2 * k))
               * 16 ** k * 5 ** (m + 1 - 2 * k) for k in range((m + 1) // 2 + 1))


def d6_model_extra(m):
    d = 6
    extra = 0
    for ct in (2, 3, 4):
        for rest in itertools.product((0, 1, 5), repeat=m):
            if (ct + sum(rest)) % d:
                continue
            c = (ct,) + rest
            sigma = sum(c) // d
            n0 = rest.count(0)
            L = (m - sigma, n0 + sigma - 1)
            base = [L[0], L[1]]
            mult = 1
            for cs in rest:
                if cs == 1:
                    base[1] += 1
                    mult *= 4
                elif cs == 5:
                    base[0] += 1
                    mult *= 4
                else:
                    mult *= 5
            for au in range(1, 6):
                for av in range(1, 6):
                    if (au + av) % d != ct:
                        continue
                    ay = (-au - av) % d
                    if ay == 0:
                        continue
                    t = (1, 0) if au + av + ay == 6 else (0, 1)
                    if (base[0] + t[0], base[1] + t[1]) == (m, m):
                        extra += mult
    return extra


def main():
    out = {"d5": [], "d6": []}
    for m in range(1, 5):
        for j in range(1, m + 2):
            n = 2 * j - 2
            Nprim = 4 if n == 0 else shioda_count(5, n)
            model = 1 + 4 ** (m + 1 - j) * Nprim
            direct = 1 + shioda_count(5, 2 * m, allowed=lambda c, j=j: all(x == 0 for x in c[j:]))
            out["d5"].append(dict(m=m, j=j, model=model, shioda=direct, ok=model == direct))
            print("d=5", out["d5"][-1]); sys.stdout.flush()
    for m in range(1, 5):
        model = 1 + R6(m) + d6_model_extra(m)
        direct = 1 + shioda_count(6, 2 * m, allowed=lambda c: all(x in (0, 1, 5) for x in c[1:]))
        out["d6"].append(dict(m=m, model=model, shioda=direct, ok=model == direct,
                              very_general=1 + R6(m), all_fermat=1 + shioda_count(6, 2 * m) if m <= 3 else None))
        print("d=6", out["d6"][-1]); sys.stdout.flush()
    with open(sys.argv[1], "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
