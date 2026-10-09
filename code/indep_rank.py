"""
indep_rank.py -- dimension counts for the split hypersurfaces (Sections 2 and 12 of the note); independent of producer.py

Object: X_b = {F_0(u_0,v_0)+...+F_m(u_m,v_m)=0} in P^{2m+1}, F_s squarefree binary forms of degree d.

Method (different from the Burau-matrix computation of producer.py):
  * Griffiths: H^{2m-q,q}_prim(X) = R_f^{(q+1)d-2m-2}, R_f = (x) R_{F_s} (Jacobian rings of the binary forms,
    Hilbert function h(e) = coeff of t^e in ((1-t^{d-1})/(1-t))^2, e in [0, 2d-4]).
  * G = mu_d^{m+1} acts by scaling (u_s,v_s); the class P*Omega/f^{q+1} with P of multidegree (e_s) has
    character c_s = e_s + 2 (mod d)  (Omega scales by zeta_s^2 per pair).  Sector of a multidegree = (c_s).
  * Single-factor piece P(c), c != 0: lower part R^{c-2} (dim c-1), upper part R^{d+c-2} (dim d-1-c).
    It is "definite" iff one part is 0 (c in {1,d-1}).  Monodromy of P(c) over the space of squarefree forms
    is FINITE iff P(u*c) is definite for every unit u mod d (definite hermitian Z[zeta]-lattice => finite
    unitary group); otherwise it is infinite and acts irreducibly (Section 3 of the
    note) -> no invariant vectors.  c = 0: P(0) = H^0_prim(d points), finite (S_d) monodromy.
  * Very general Hodge rank R(d,m) = sum over multidegrees of total degree (m+1)(d-2) whose sector lies in
    F_d^{m+1} of prod h(e_s).
Checks:
  (i)  sum of all pieces = b_prim = ((d-1)^{2m+2} + (d-1))/d, and the (m,m) part = h^{m,m}_prim computed
       independently by counting Fermat characters (a_0..a_{2m+1}), a_i in [1,d-1], sum = (m+1)d.
  (ii) R(d,m) vs. closed formula sum_k (m+1)!/(k!^2 (m+1-2k)!) (d-2)^{2k} (d-1)^{m+1-2k} (d in {3,4,6})
       and (d-1)^{m+1} otherwise; vs. values from other programs and Lin-Zhang's (d-1)^{m+1}.
  (iii) every F_d-sector piece that meets H^{m,m} is entirely of type (m,m) iff balanced (#1 = #(d-1));
       non-balanced F_d-sectors contain no (m,m) part.
  (iv) first-order IVHS upper bound (classes killed by all split deformations) >= R (sanity only).
  (v)  surface case m=1: rho = 1 + R(d,1) vs. cited classical values (Sasakura 1968 prime d: (d-1)^2+1;
       Shioda 1981 eq. (5.11)/Rem. 5.3: generic quartic P(x,y)+Q(z,w): 18 ... ).
"""
import itertools, json, math
from math import factorial


def hilb_binary(d):
    # coefficients of ((1 - t^(d-1)) / (1 - t))^2 = (1 + t + ... + t^(d-2))^2
    h = [0] * (2 * d - 3)
    for i in range(d - 1):
        for j in range(d - 1):
            h[i + j] += 1
    return h  # index e = 0 .. 2d-4


def finite_set(d):
    F = {0}
    for c in range(1, d):
        orbit = {(u * c) % d for u in range(1, d) if math.gcd(u, d) == 1}
        definite = all(oc in (1, d - 1) for oc in orbit)
        if definite:
            F.add(c)
    return F


def closed_formula(d, m):
    if d in (3, 4, 6):
        return sum(factorial(m + 1) // (factorial(k) ** 2 * factorial(m + 1 - 2 * k))
                   * (d - 2) ** (2 * k) * (d - 1) ** (m + 1 - 2 * k) for k in range(0, (m + 1) // 2 + 1))
    return (d - 1) ** (m + 1)


def fermat_hmm(d, m):
    # independent count of h^{m,m}_prim: Fermat characters with sum (m+1)d  (Griffiths/Shioda)
    n2 = 2 * m + 2
    # dynamic programming over number of variables
    dp = {0: 1}
    for _ in range(n2):
        nd = {}
        for s, c in dp.items():
            for a in range(1, d):
                nd[s + a] = nd.get(s + a, 0) + c
        dp = nd
    return dp.get((m + 1) * d, 0)


def analyse(d, m):
    h = hilb_binary(d)
    F = finite_set(d)
    E_mm = (m + 1) * (d - 2)
    tot = 0
    hmm = 0
    R = 0
    sector_report = {}
    nonbal_mm = 0
    for es in itertools.product(range(2 * d - 3), repeat=m + 1):
        E = sum(es)
        if (E + 2 * m + 2) % d:
            continue
        dim = 1
        for e in es:
            dim *= h[e]
        if dim == 0:
            continue
        tot += dim
        c = tuple((e + 2) % d for e in es)
        assert sum(c) % d == 0
        if E == E_mm:
            hmm += dim
            if all(x in F for x in c):
                R += dim
                n1 = sum(1 for x in c if x == 1 and d > 2 and 1 != d - 1)
                nm = sum(1 for x in c if x == d - 1 and 1 != d - 1)
                if n1 != nm:
                    nonbal_mm += dim
        if all(x in F for x in c):
            key = 'bal' if sum(1 for x in c if x == 1) == sum(1 for x in c if x == d - 1) else 'nonbal'
            q = (E + 2 * m + 2) // d - 1
            sector_report.setdefault(key, {}).setdefault(str((2 * m - q, q)), 0)
            sector_report[key][str((2 * m - q, q))] += dim
    bprim = ((d - 1) ** (2 * m + 2) + (d - 1)) // d
    # IVHS first-order kernel: multidegrees of degree E_mm with every e_s >= d-3 (Ann_s(e)=R^e iff e >= d-3)
    ivhs = 0
    for es in itertools.product(range(2 * d - 3), repeat=m + 1):
        if sum(es) == E_mm and all(e >= d - 3 for e in es):
            dim = 1
            for e in es:
                dim *= h[e]
            ivhs += dim
    return dict(d=d, m=m, finite_sectors=sorted(F), b_prim_sum=tot, b_prim_formula=bprim,
                hmm_prim_jacobian=hmm, hmm_prim_fermat=fermat_hmm(d, m), R=R,
                R_closed=closed_formula(d, m), nonbalanced_in_mm=nonbal_mm,
                finite_sector_hodge_types=sector_report, ivhs_first_order_bound=ivhs)


if __name__ == '__main__':
    out = {}
    ok = True
    for d in range(3, 10):
        for m in range(1, 5):
            if (2 * d - 3) ** (m + 1) > 3_000_000:
                continue
            r = analyse(d, m)
            out[f'{d},{m}'] = r
            ok &= (r['b_prim_sum'] == r['b_prim_formula'])
            ok &= (r['hmm_prim_jacobian'] == r['hmm_prim_fermat'])
            ok &= (r['R'] == r['R_closed'])
            ok &= (r['nonbalanced_in_mm'] == 0)
            ok &= (r['ivhs_first_order_bound'] >= r['R'])
            print(f"d={d} m={m}  F_d={r['finite_sectors']}  R={r['R']} (closed {r['R_closed']})  "
                  f"h^mm_prim={r['hmm_prim_jacobian']}/{r['hmm_prim_fermat']}  b_prim={r['b_prim_sum']}/{r['b_prim_formula']}  "
                  f"IVHS1>= {r['ivhs_first_order_bound']}  types(finite sectors)={r['finite_sector_hodge_types']}")
    # comparison with values from other programs (producer.py, and one not included here) and Lin-Zhang's (d-1)^{m+1}
    claimed = {'3,2': 20, '3,3': 70, '4,2': 99, '4,3': 609, '4,4': 3843, '6,2': 605, '6,3': 6961, '6,4': 81525,
               '5,2': 64, '5,3': 256, '7,2': 216, '8,2': 343, '7,3': 1296, '8,3': 2401, '9,4': 32768,
               '3,1': 6, '4,1': 17, '5,1': 16, '6,1': 57, '7,1': 36, '8,1': 49}
    cmp = {}
    for k, v in claimed.items():
        d, m = map(int, k.split(','))
        mine = out[k]['R'] if k in out else closed_formula(d, m)
        src = 'enumerated' if k in out else 'closed formula only (enumeration too large)'
        cmp[k] = dict(claimed=v, mine=mine, agree=(mine == v), how=src)
        ok &= (mine == v)
    # surface (m=1) classical controls: rho = 1 + R(d,1)
    surf = {d: 1 + out[f'{d},1']['R'] for d in range(3, 10)}
    classical = {3: (7, 'every smooth cubic surface: rho=7'),
                 4: (18, 'generic P(x,y)+Q(z,w) quartic K3: Shioda 1981 (5.11) with r=8; cf. Rem 5.3(ii): P=Q gives 19'),
                 5: (17, 'Sasakura 1968 (prime m): (m-1)^2+1 for generic P,Q (cited via Shioda 1981 Rem 5.3(i))'),
                 7: (37, 'Sasakura 1968 (prime m): (m-1)^2+1')}
    for d, (v, why) in classical.items():
        ok &= (surf[d] == v)
    out['_comparison_with_claims'] = cmp
    out['_surface_rho'] = {str(d): dict(mine=surf[d], classical=classical.get(d)) for d in surf}
    out['_all_checks_pass'] = ok
    print('comparison:', json.dumps(cmp))
    print('surface rho:', surf)
    print('ALL CHECKS PASS:', ok)
    with open('indep_rank_out.json', 'w') as fh:
        json.dump(out, fh, indent=1)
