"""Control for i3_lines.py: a sublattice of index 2 (one split line doubled) must give index^2 = 4."""
import sys
sys.path.insert(0, sys.argv[1])
import i3_lines as L
M = L.gram(); n = len(L.LINES)
split = [i for i, (p, a, b) in enumerate(L.LINES) if p == 0]
rows = [M[i] for i in split]
B = L.hnf_rows(rows)
B2 = [list(B[0])] + [list(r) for r in B[1:]]
B2[0] = [2 * x for x in B2[0]]
bl = []
for i in range(n):
    if L.rank_q([M[j] for j in bl + [i]]) > len(bl):
        bl.append(i)
    if len(bl) == 37: break
from fractions import Fraction
Gb = [[Fraction(M[i][j]) for j in bl] for i in bl]
def coords(v):
    A = [[Gb[j][i] for j in range(37)] + [Fraction(v[bl[i]])] for i in range(37)]
    for c in range(37):
        p = next(i for i in range(c, 37) if A[i][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for i in range(37):
            if i != c and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [A[i][37] for i in range(37)]
def gram(rows):
    C = [coords(v) for v in rows]
    return [[sum(C[x][i] * Gb[i][j] * C[y][j] for i in range(37) for j in range(37)) for y in range(len(rows))] for x in range(len(rows))]
d1 = L.det_int(gram(B)); d2 = L.det_int(gram(B2))
print("disc split", d1, "disc doubled", d2, "ratio", d2 / d1, "(expect 4)")
