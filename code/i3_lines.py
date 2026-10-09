"""i3_lines.py -- integral Neron-Severi of the very general split quintic surface, checked at the Fermat member.

X = {u0^5+v0^5+u1^5+v1^5 = 0} = {x0^5+x1^5+x2^5+x3^5=0}, pairs (x0,x1),(x2,x3).
75 lines: for each of the 3 pairings {i,j}|{k,l}: x_i = a x_j, x_k = b x_l with a^5 = b^5 = -1.
Split lines (pairing 01|23) are the Lambda_p of the note (m=1), 25 of them.
Intersection numbers: L.L = 2-d = -3; L.L' = 1 if the lines meet, 0 otherwise (decided exactly in Q(zeta_10)).
G = mu_5^2 acts by (x0,x1,x2,x3) -> (g0 x0, g0 x1, g1 x2, g1 x3); it permutes the lines.

Facts used (not recomputed): Degtyarev, Thm 1.2: NS(Fermat quintic surface) = Z-span of the 75 lines.
For very general F (d=5, m=1): NS(X_F) = H^2(X_F,Z) cap V_0 = H^2(X_F,Z)^G, a locally constant lattice;
at the Fermat member it equals NS(Fermat)^G = (Z-span of lines)^G.
Checks:
  rank(span lines) = 37; rank(span split lines) = 17;
  N^G (G-invariants of N = Z-span of all lines, computed as an abstract lattice Z^75/radical);
  index [N^G : span_Z(split lines)] (expected 1 by Degtyarev-Shimada Cor. 1.6);
  discriminants. Plain Python, exact integers/Fractions.
"""
import itertools, json, sys
from fractions import Fraction

# ---- Q(w), w = primitive 10th root of unity, Phi_10 = x^4 - x^3 + x^2 - x + 1; elements: tuples of 4 Fractions/ints
def mul(a, b):
    r = [0] * 7
    for i in range(4):
        if a[i] == 0:
            continue
        for j in range(4):
            r[i + j] += a[i] * b[j]
    # reduce x^k, k>=4, using x^4 = x^3 - x^2 + x - 1
    for k in range(6, 3, -1):
        c = r[k]
        if c:
            r[k] = 0
            r[k - 1] += c
            r[k - 2] -= c
            r[k - 3] += c
            r[k - 4] -= c
    return tuple(r[:4])


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)


def wpow(k):
    r = ONE
    for _ in range(k % 10):
        r = mul(r, (0, 1, 0, 0))
    return r


def det(M):
    """Determinant of a square matrix over Z[w] by cofactor expansion (size 4)."""
    n = len(M)
    if n == 1:
        return M[0][0]
    tot = ZERO
    for j in range(n):
        if M[0][j] == ZERO:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        t = mul(M[0][j], det(minor))
        tot = add(tot, t if j % 2 == 0 else neg(t))
    return tot


ROOTS = [wpow(k) for k in (1, 3, 5, 7, 9)]  # a^5 = -1  (w^5 = -1)
assert all(mul(mul(mul(mul(a, a), a), a), a) == neg(ONE) for a in ROOTS)
PAIRINGS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]


def line_basis(pairing, ia, ib):
    """Two spanning vectors of {x_i = a x_j, x_k = b x_l}: (x_j=1 => x_i=a), (x_l=1 => x_k=b)."""
    (i, j), (k, l) = pairing
    a, b = ROOTS[ia], ROOTS[ib]
    v1 = [ZERO] * 4
    v1[j] = ONE
    v1[i] = a
    v2 = [ZERO] * 4
    v2[l] = ONE
    v2[k] = b
    return v1, v2


LINES = [(pi, ia, ib) for pi in range(3) for ia in range(5) for ib in range(5)]


def on_fermat(v):
    s = ZERO
    for x in v:
        x5 = mul(mul(mul(mul(x, x), x), x), x)
        s = add(s, x5)
    return s == ZERO


def gram():
    vecs = [line_basis(PAIRINGS[p], a, b) for (p, a, b) in LINES]
    for (v1, v2) in vecs:  # the lines lie on X: check points v1, v2, v1+v2, v1+w v2
        for v in (v1, v2, [add(x, y) for x, y in zip(v1, v2)], [add(x, mul((0, 1, 0, 0), y)) for x, y in zip(v1, v2)]):
            assert on_fermat(v)
    n = len(LINES)
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = -3
        for j in range(i + 1, n):
            A = [vecs[i][0], vecs[i][1], vecs[j][0], vecs[j][1]]
            meet = det(A) == ZERO
            M[i][j] = M[j][i] = 1 if meet else 0
    return M


def g_action(g0, g1):
    """Permutation of LINES induced by (x0,x1,x2,x3) -> (z^g0 x0, z^g0 x1, z^g1 x2, z^g1 x3), z = w^2 (5th root)."""
    perm = []
    scale = {0: g0, 1: g0, 2: g1, 3: g1}
    for (p, ia, ib) in LINES:
        (i, j), (k, l) = PAIRINGS[p]
        # image line: x_i' = a * z^{s_i - s_j} x_j'
        ka = (2 * ia + 1 + 2 * (scale[i] - scale[j])) % 10
        kb = (2 * ib + 1 + 2 * (scale[k] - scale[l])) % 10
        perm.append(LINES.index((p, (ka - 1) // 2, (kb - 1) // 2)))
    return perm


# ---- integer linear algebra
def hnf_rows(rows):
    """Row-style Hermite normal form (echelon over Z) of a list of integer rows; returns nonzero rows."""
    A = [list(r) for r in rows]
    m = len(A)
    ncol = len(A[0]) if A else 0
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, m):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        while True:
            # Euclid on column c among rows r..m-1
            nz = [i for i in range(r + 1, m) if A[i][c] != 0]
            if not nz:
                break
            for i in nz:
                q = A[i][c] // A[r][c]
                A[i] = [x - q * y for x, y in zip(A[i], A[r])]
            nz = [i for i in range(r, m) if A[i][c] != 0]
            best = min(nz, key=lambda i: abs(A[i][c]))
            A[r], A[best] = A[best], A[r]
        if A[r][c] < 0:
            A[r] = [-x for x in A[r]]
        for i in range(r):
            q = A[i][c] // A[r][c]
            A[i] = [x - q * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == m:
            break
    return [row for row in A[:r]]


def det_int(M):
    M = [[Fraction(x) for x in row] for row in M]
    n = len(M)
    d = Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            if f:
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def rank_q(M):
    M = [[Fraction(x) for x in row] for row in M]
    r = 0
    ncol = len(M[0])
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r


def main():
    M = gram()
    n = len(LINES)
    split = [i for i, (p, a, b) in enumerate(LINES) if p == 0]
    rk_all = rank_q(M)
    rk_split = rank_q([[M[i][j] for j in range(n)] for i in split])
    # Abstract lattice N = Z^75 / radical. Coordinates of a line class: its row of intersection numbers with all lines
    # (the map N -> Z^75, x -> (x.L_j) is injective because the form is nondegenerate on N).
    # Z-basis of the image: HNF of the 75 rows. Elements of N <-> integer combinations of rows.
    B = hnf_rows(M)  # rows: Z-basis of the image of N (rank 37)
    assert len(B) == rk_all
    # G-invariants of N: x in N invariant iff its row vector v satisfies v[perm(j)] = v[j]  (since g preserves the form
    # and permutes lines: (g x).L_{g j} = x.L_j). Solve over Z: v = sum_i c_i B_i, constraints linear in c.
    gens = [g_action(1, 0), g_action(0, 1)]
    for P in gens:  # check G preserves the Gram matrix
        assert all(M[P[i]][P[j]] == M[i][j] for i in range(n) for j in range(n))
    for i in split:
        assert all(P[i] == i for P in gens)
    # constraint rows: for each generator P and each j: v[P[j]] - v[j] = 0  ->  sum_i c_i (B_i[P[j]] - B_i[j]) = 0
    k = len(B)
    cons = []
    for P in gens:
        for j in range(n):
            cons.append([B[i][P[j]] - B[i][j] for i in range(k)])
    # integer kernel of cons (as c in Z^k): use HNF on [cons^T | I]
    aug = [[cons[r][i] for r in range(len(cons))] + [1 if t == i else 0 for t in range(k)] for i in range(k)]
    H = hnf_rows(aug)
    ncons = len(cons)
    # rows of the HNF whose constraint part is zero give a Z-basis of the kernel (the transform is unimodular)
    full = [list(r) for r in H]
    # hnf_rows drops zero rows only; rows with zero constraint part are kept since identity part nonzero
    ker = [r[ncons:] for r in full if all(x == 0 for x in r[:ncons])]
    inv_rows = [[sum(c[i] * B[i][j] for i in range(k)) for j in range(n)] for c in ker]
    rk_inv = rank_q(inv_rows)
    # Gram on invariant lattice: for x,y in N given by row vectors v_x, v_y, x.y = ? We need x as a combination of
    # lines: solve. Easier: express everything in the coordinates of a Q-basis of 37 lines.
    # choose 37 independent lines
    basis_lines = []
    for i in range(n):
        if rank_q([M[j] for j in basis_lines + [i]]) > len(basis_lines):
            basis_lines.append(i)
        if len(basis_lines) == rk_all:
            break
    Gb = [[Fraction(M[i][j]) for j in basis_lines] for i in basis_lines]
    # coordinates of a vector with intersection row v: solve Gb^T c = (v restricted to basis_lines)
    def coords(v):
        A = [row[:] + [Fraction(v[basis_lines[t]])] for t, row in enumerate([[Gb[j][i] for j in range(rk_all)] for i in range(rk_all)])]
        # Gauss-Jordan
        m_ = rk_all
        for c in range(m_):
            p = next(i for i in range(c, m_) if A[i][c] != 0)
            A[c], A[p] = A[p], A[c]
            pv = A[c][c]
            A[c] = [x / pv for x in A[c]]
            for i in range(m_):
                if i != c and A[i][c] != 0:
                    f = A[i][c]
                    A[i] = [x - f * y for x, y in zip(A[i], A[c])]
        return [A[i][m_] for i in range(m_)]
    def gram_of(rows):
        C = [coords(v) for v in rows]
        return [[sum(a * b * 0 for a, b in []) + sum(C[x][i] * Gb[i][j] * C[y][j] for i in range(rk_all) for j in range(rk_all))
                 for y in range(len(rows))] for x in range(len(rows))]
    G_inv = gram_of(inv_rows)
    disc_inv = det_int(G_inv)
    # split lattice: rows of the split lines; Z-basis via HNF
    S_rows = hnf_rows([M[i] for i in split])
    G_split = gram_of(S_rows)
    disc_split = det_int(G_split)
    disc_all = det_int(gram_of(B))
    idx2 = disc_split / disc_inv
    out = dict(rank_all_lines=rk_all, rank_split=rk_split, rank_G_invariants=rk_inv,
               disc_NS_fermat=str(disc_all), disc_NS_G_invariant=str(disc_inv), disc_split_span=str(disc_split),
               index_squared=str(idx2))
    print(json.dumps(out, indent=1))
    with open(sys.argv[1] if len(sys.argv) > 1 else "i3_lines.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
