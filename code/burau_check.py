"""burau_check.py: second implementation of the monodromy checks, plain Python, own arithmetic in Z[zeta_d].
Model (different from burau_sage.py): the reduced Burau representation in the convention of Lin-Zhang (arXiv 2609.03784,
Sec. 6.3), specialised at t = zeta_d^r, on V = Q(zeta_d)^(d-1):
   sigma_i: e_i -> -t e_i,  e_(i-1) -> e_(i-1) + t e_i,  e_(i+1) -> e_(i+1) + e_i,  e_j -> e_j otherwise,
the invariant vector v_t = sum_j (1 + t + ... + t^(j-1)) e_j and the quotient E = V / C v_t (dim d - 2), with the
hermitian form H_t of (6.6) (r != d/2) or (zeta - zeta^-1) times the alternating form (6.11) (r = d/2);
the form is h(x, y) = x^T H conj(y).
For d = 4, 5, 6 and every r: invariance of the form, B_i v_t = v_t, braid relations, signature on E, dimension of the
algebra generated on E, an infinite-order witness (|sigma(trace)| > d - 2 for some Galois embedding sigma) when some
conjugate is indefinite, and the order of the image group (breadth-first closure) when every conjugate is definite.
usage: python burau_check.py
"""
import json, cmath, math, itertools, time
from fractions import Fraction as Fr

PHI = {4: [1, 0, 1], 5: [1, 1, 1, 1, 1], 6: [1, -1, 1], 3: [1, 1, 1]}     # coefficients, low degree first, monic


class Cyc:
    """Q(zeta_d) as polynomials mod Phi_d (coefficients Fraction or int)"""
    def __init__(self, d):
        self.d, self.phi = d, PHI[d]
        self.n = len(self.phi) - 1

    def red(self, a):
        a = list(a) + [0] * max(0, self.n - len(a))
        for k in range(len(a) - 1, self.n - 1, -1):
            c = a[k]
            if c:
                for j in range(self.n + 1):
                    a[k - self.n + j] -= c * self.phi[j]
        return tuple(a[:self.n])

    def zeta(self, e):
        a = [0] * (e % self.d + 1)
        a[e % self.d] = 1
        return self.red(a)

    def const(self, c):
        return self.red([c])

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def mul(self, a, b):
        r = [0] * (2 * self.n - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        r[i + j] += x * y
        return self.red(r)

    def conj(self, a):
        r = [0] * self.d
        for i, x in enumerate(a):
            r[(-i) % self.d] += x
        return self.red(r)

    def iszero(self, a):
        return all(x == 0 for x in a)

    def inv(self, a):
        # extended Euclid in Q[x]
        def trim(p):
            p = list(p)
            while p and p[-1] == 0:
                p.pop()
            return p

        def pdivmod(p, q):
            p = [Fr(x) for x in p]
            q = trim(q)
            out = [Fr(0)] * max(1, len(p) - len(q) + 1)
            while len(trim(p)) >= len(q):
                p = trim(p)
                c = p[-1] / q[-1]
                s = len(p) - len(q)
                out[s] = c
                for j in range(len(q)):
                    p[s + j] -= c * q[j]
            return trim(out), trim(p)

        def pmul(p, q):
            r = [Fr(0)] * (len(p) + len(q) - 1) if p and q else []
            for i, x in enumerate(p):
                for j, y in enumerate(q):
                    r[i + j] += x * y
            return r

        def psub(p, q):
            n = max(len(p), len(q))
            return trim([(p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0) for i in range(n)])
        r0, r1 = trim([Fr(x) for x in self.phi]), trim([Fr(x) for x in a])
        s0, s1 = [], [Fr(1)]
        while r1:
            qq, rr = pdivmod(r0, r1)
            r0, r1 = r1, rr
            s0, s1 = s1, psub(s0, pmul(qq, s1))
        assert len(r0) == 1, 'not invertible'
        return self.red([x / r0[0] for x in s0] + [0] * self.n)

    def emb(self, a, u):
        w = cmath.exp(2j * math.pi * u / self.d)
        return sum(complex(x) * w ** i for i, x in enumerate(a))


def mat_mul(F, A, B):
    n, m, p = len(A), len(B), len(B[0])
    Z = F.const(0)
    out = []
    for i in range(n):
        row = []
        for j in range(p):
            s = Z
            for k in range(m):
                if not F.iszero(A[i][k]) and not F.iszero(B[k][j]):
                    s = F.add(s, F.mul(A[i][k], B[k][j]))
            row.append(s)
        out.append(row)
    return out


def ident(F, n):
    return [[F.const(1 if i == j else 0) for j in range(n)] for i in range(n)]


def conjT(F, A):
    return [[F.conj(A[j][i]) for j in range(len(A))] for i in range(len(A[0]))]


def rank(F, rows):
    A = [list(r) for r in rows]
    r = 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if not F.iszero(A[i][c])), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = F.inv(A[r][c])
        A[r] = [F.mul(x, inv) for x in A[r]]
        for i in range(len(A)):
            if i != r and not F.iszero(A[i][c]):
                f = A[i][c]
                A[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(A[i], A[r])]
        r += 1
    return r


def burau(F, d, r):
    t = F.zeta(r)
    n = d - 1
    Bs = []
    for i in range(n):
        B = ident(F, n)
        B[i][i] = F.neg(t)
        if i > 0:
            B[i][i - 1] = t
        if i < n - 1:
            B[i][i + 1] = F.const(1)
        Bs.append(B)
    return t, Bs


def form(F, d, r, t):
    n = d - 1
    H = [[F.const(0)] * n for _ in range(n)]
    if 2 * r != d:
        tinv = F.inv(t)
        for j in range(n):
            H[j][j] = F.add(F.const(2), F.add(t, tinv))
            if j + 1 < n:
                H[j][j + 1] = F.neg(F.add(F.const(1), t))
                H[j + 1][j] = F.neg(F.add(F.const(1), tinv))
    else:
        s = F.sub(F.zeta(1), F.zeta(-1))                            # zeta - zeta^-1: purely imaginary, nonzero
        for j in range(n - 1):
            H[j][j + 1] = F.mul(s, F.const(-1))
            H[j + 1][j] = s
    return H


def quotient_mats(F, Bs, v):
    """E = V / C v with complement e_2..e_n (v has e_1-coefficient 1): x -> x - x_1 v"""
    n = len(v)
    assert v[0] == F.const(1)
    out = []
    for B in Bs:
        A = []
        for i in range(1, n):
            row = []
            for j in range(1, n):
                # image of e_j (column j of B), projected: coefficient i of (Be_j - (Be_j)_1 v)
                c = F.sub(B[i][j], F.mul(B[0][j], v[i]))
                row.append(c)
            A.append(row)
        out.append(A)
    return out


def herm_signature(F, H):
    H = [list(r) for r in H]
    idx = list(range(len(H)))
    pos = neg = 0
    while idx:
        piv = next((i for i in idx if not F.iszero(H[i][i])), None)
        if piv is None:
            pair = next(((i, j) for i in idx for j in idx if i < j and not F.iszero(H[i][j])), None)
            if pair is None:
                break
            i, j = pair
            for c in (F.const(1), F.zeta(1), F.zeta(2)):
                val = F.add(F.mul(c, H[j][i]), F.mul(F.conj(c), H[i][j]))
                if not F.iszero(val):
                    break
            # e_i -> e_i + c e_j : row i += c * row j ; column i += conj(c) * column j
            H[i] = [F.add(x, F.mul(c, y)) for x, y in zip(H[i], H[j])]
            for k in range(len(H)):
                H[k][i] = F.add(H[k][i], F.mul(F.conj(c), H[k][j]))
            continue
        p = H[piv][piv]
        val = F.emb(p, 1)
        assert abs(val.imag) < 1e-9 and abs(val.real) > 1e-9
        if val.real > 0:
            pos += 1
        else:
            neg += 1
        pinv = F.inv(p)
        for i in idx:
            if i != piv and not F.iszero(H[i][piv]):
                f = F.mul(H[i][piv], pinv)
                for j in idx:
                    H[i][j] = F.sub(H[i][j], F.mul(f, H[piv][j]))
        idx.remove(piv)
    return pos, neg


def algebra_dim(F, As):
    m = len(As[0])
    span = []
    I = ident(F, m)

    def add(M):
        v = [x for row in M for x in row]
        if rank(F, span + [v]) > len(span):
            span.append(v)
            return True
        return False
    add(I)
    frontier = [I]
    while frontier and len(span) < m * m:
        new = []
        for M in frontier:
            for A in As:
                P = mat_mul(F, M, A)
                if add(P):
                    new.append(P)
        frontier = new
    return len(span)


def inverse(F, A):
    n = len(A)
    M = [list(A[i]) + [F.const(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(i for i in range(c, n) if not F.iszero(M[i][c]))
        M[c], M[piv] = M[piv], M[c]
        inv = F.inv(M[c][c])
        M[c] = [F.mul(x, inv) for x in M[c]]
        for i in range(n):
            if i != c and not F.iszero(M[i][c]):
                f = M[i][c]
                M[i] = [F.sub(x, F.mul(f, y)) for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]


def witness(F, As, units, maxlen=5):
    gens = [('s%d' % (i + 1), A) for i, A in enumerate(As)] + [('s%d^-1' % (i + 1), inverse(F, A)) for i, A in enumerate(As)]
    m = len(As[0])
    for L in range(1, maxlen + 1):
        for word in itertools.product(range(len(gens)), repeat=L):
            M = ident(F, m)
            for w in word:
                M = mat_mul(F, M, gens[w][1])
            tr = F.const(0)
            for i in range(m):
                tr = F.add(tr, M[i][i])
            for u in units:
                if abs(F.emb(tr, u)) > m + 1e-6:
                    return dict(word='*'.join(gens[w][0] for w in word), embedding_unit=u, abs_trace=abs(F.emb(tr, u)))
    return None


def group_order(F, As, cap=400000):
    key = lambda M: tuple(int(x) for row in M for e in row for x in e)
    I = ident(F, len(As[0]))
    seen = {key(I)}
    frontier = [I]
    while frontier:
        new = []
        for M in frontier:
            for A in As:
                P = mat_mul(F, M, A)
                k = key(P)
                if k not in seen:
                    seen.add(k)
                    new.append(P)
                    if len(seen) > cap:
                        return None
        frontier = new
    return len(seen)


def main():
    out = []
    for d in (4, 5, 6):
        F = Cyc(d)
        units = [u for u in range(1, d) if math.gcd(u, d) == 1]
        for r in range(1, d):
            t0 = time.time()
            t, Bs = burau(F, d, r)
            n = d - 1
            H = form(F, d, r, t)
            # the form is h(x, y) = x^T H conj(y); invariance: B^T H conj(B) = H  (Lin-Zhang (6.7) write B^* H B = H,
            # which holds for H^T = conj(H) instead: same signature and radical)
            BT = lambda B: [[B[j][i] for j in range(len(B))] for i in range(len(B))]
            Bc = lambda B: [[F.conj(x) for x in row] for row in B]
            inv_ok = all(mat_mul(F, mat_mul(F, BT(B), H), Bc(B)) == H for B in Bs)
            v = []
            acc = F.const(0)
            for j in range(n):
                acc = F.add(acc, F.zeta(r * j))
                v.append(acc)
            fix_ok = all([row[0] for row in mat_mul(F, B, [[x] for x in v])] == v for B in Bs)
            braid_ok = all(mat_mul(F, mat_mul(F, Bs[i], Bs[i + 1]), Bs[i]) == mat_mul(F, mat_mul(F, Bs[i + 1], Bs[i]), Bs[i + 1])
                           for i in range(n - 1)) and \
                all(mat_mul(F, Bs[i], Bs[j]) == mat_mul(F, Bs[j], Bs[i]) for i in range(n) for j in range(n) if abs(i - j) > 1)
            As = quotient_mats(F, Bs, v)
            HE = [row[1:] for row in H[1:]]                                   # v in the radical: form on e_2..e_n
            rad_ok = all(F.iszero(F.add(F.const(0), sum_row)) for sum_row in
                         [[x for x in mat_mul(F, [v], H)[0]][k] for k in range(n)])
            sig = herm_signature(F, HE)
            conj = sorted(set((u * r) % d for u in units))
            indef = [c for c in conj if 2 <= c <= d - 2]
            alg = algebra_dim(F, As)
            res = dict(d=d, r=r, form_invariant=inv_ok, v_fixed=fix_ok, braid_relations=braid_ok, v_in_radical=rad_ok,
                       dim_E=len(As[0]), signature=sig, algebra_dim=alg, absolutely_irreducible=alg == len(As[0]) ** 2,
                       galois_conjugates=conj, indefinite_conjugates=indef)
            if indef:
                res['infinite_order_witness'] = witness(F, As, units)
            else:
                res['group_order'] = group_order(F, As)
            res['seconds'] = round(time.time() - t0, 1)
            out.append(res)
            print(json.dumps(res), flush=True)
    json.dump(out, open('burau_check.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
