Supporting code and data for the note
"The Hodge conjecture for very general sums of binary quintics and sextics"
==========================================================================

Status: draft accompanying a draft note; not peer reviewed.

What the note proves (Theorem A there): for d = 5, 6 and every m >= 1, the Hodge conjecture holds for
X = {F_0(u_0,v_0) + ... + F_m(u_m,v_m) = 0} in P^(2m+1) when F_0, ..., F_m are very general binary forms
of degree d, and dim Hdg^m(X) = 1 + R(d,m) with R(5,m) = 4^(m+1) and
R(6,m) = sum_k (m+1)!/(k!^2 (m+1-2k)!) 16^k 5^(m+1-2k).

IMPORTANT: the proofs in the note do not use any of these computations. The scripts below are
consistency checks of the dimension counts (Section 2 of the note) and of the monodromy statements
(Section 3), for finitely many cases; they are described in Section 8 of the note. The numbers quoted
in the note (values of R(d,m), h^{2,2} = 1752 for d = 6, m = 2) are evaluations of formulas stated
there; indep_rank.py recomputes them.

Licence: the code in code/ is under the MIT licence; the data in outputs/ and this text are under
CC BY 4.0 (see the LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3; numpy only for producer.py, the other Python scripts use only the standard library.
burau_sage.py needs SageMath (it uses GAP, included in SageMath, for the group orders). Tested with
Python 3.12.9 and numpy 2.1.3, and with SageMath 10.10 in the docker image sagemath/sagemath
(digest sha256:310dfe23fc786f3a5a06a7b0387c60201ec79b32c0db941b7624505ee5edb565), limited to 2 CPUs
and 3 GB of memory. The scripts do not parallelize themselves; numpy may start several threads, which
can be prevented by setting OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and MKL_NUM_THREADS=1 (this was
done for the runtimes below). Runtimes measured on a desktop PC, one script at a time (another
computation was running on the machine): indep_rank.py 1 s, producer.py 71 s, check.py 33 s,
burau_check.py 104 s, burau_sage.py 14 s (docker, including start-up).


2. How to run
-------------
Each script writes its result file (.json) into the CURRENT working directory and prints a summary
to the screen. The reference results are in outputs/ (the .out files there are the screen output).
To reproduce them, from the top folder of this repository:

    mkdir rerun
    cd rerun
    python ../code/indep_rank.py   > indep_rank.out
    python ../code/producer.py     > producer.out
    python ../code/check.py        > check.out       (reads producer_out.json, so run it in the
                                                      same folder, after producer.py)
    python ../code/burau_check.py  > burau_check.out
    sage -python ../code/burau_sage.py > burau_sage.out
    python ../code/compare_outputs.py ../outputs .

With docker instead of a local SageMath, the fifth command is, from the rerun folder (bash; on
Windows with Git Bash, prefix the docker line with MSYS_NO_PATHCONV=1):

    IMG=sagemath/sagemath@sha256:310dfe23fc786f3a5a06a7b0387c60201ec79b32c0db941b7624505ee5edb565
    docker run --rm -v "$(pwd)/..:/w" -w /w/rerun $IMG sage -python ../code/burau_sage.py > burau_sage.out

compare_outputs.py compares every .json file in outputs/ with the file of the same name in the
current folder; it ignores the run times (keys "seconds") and compares floating-point numbers with
relative tolerance 1e-9. Expected last line: "ALL EQUAL". The .out files of burau_sage.py and
burau_check.py contain run times and differ from run to run in those.

Reproducibility: the files in outputs/ were produced by the commands above in a fresh folder (Windows,
Python and numpy versions as above, SageMath in the docker image above). They agree with the results
of earlier runs of the same code, which has since changed only in comments (compare_outputs.py:
"ALL EQUAL"; four of the five .json files byte for byte, burau_check.json up to run times). All
content is exact integer, rational or cyclotomic arithmetic or uses fixed random seeds, except the
floating-point quantities of producer.py (residuals "fix_err", "inv_err", eigenvalue moduli
"max_abs_eig"), which may differ in the last digits with other numpy/BLAS versions.


3. Files
--------
Notation as in the note: d = degree of the binary forms, m + 1 = number of blocks, so X has dimension
2m; R = R_0 (x) ... (x) R_m is the Jacobian ring, each factor with Hilbert function
(1 + t + ... + t^(d-2))^2; a multidegree (e_0, ..., e_m) has sector c_s = e_s + 2 (mod d).

code/indep_rank.py  ->  outputs/indep_rank_out.json, outputs/indep_rank.out
  Supports: Section 2 of the note (Lemma 2.2, Definition 2.3 and the closed formula for R(d,m)) and
  the values quoted in Theorem A, Section 8 and Remark 7.4.
  Checks, for 3 <= d <= 9 and 1 <= m <= 4, by enumerating the multidegrees of R: (i) the dimensions
  of all pieces add up to the primitive Betti number ((d-1)^(2m+2) + d - 1)/d, and the dimension of
  the (m,m)-part agrees with an independent count of Fermat characters; (ii) R(d,m), the dimension of
  the (m,m)-part of the sectors with all c_s in Phi_d (Phi_d = {0, 1, -1} for d = 3, 4, 6 and {0}
  otherwise), equals the closed formula, and agrees with values from other programs (producer.py
  below, and one not included here) and with Lin-Zhang's (d-1)^(m+1) for d >= 7; (iii) every
  Phi_d-sector that meets H^{m,m} is balanced (as many c_s = 1 as c_s = -1) and entirely of type
  (m,m), and no other Phi_d-sector meets H^{m,m}; (iv) an upper bound from first-order deformations
  (sanity check only, not used); (v) m = 1: rho = 1 + R(d,1) equals the classical Picard numbers 7
  (d = 3), 18 (d = 4), 17 (d = 5), 37 (d = 7).
  Expected: "ALL CHECKS PASS: True" and "_all_checks_pass": true; for example R(5,m) = 16, 64, 256,
  1024 and R(6,m) = 57, 605, 6961, 81525 for m = 1, 2, 3, 4; for (d,m) = (6,2) the (m,m)-part has
  dimension 1751 (h^{2,2} = 1752 with the square of the hyperplane class).

code/producer.py  ->  outputs/producer_out.json, outputs/producer.out
  Supports: Section 8 (a second computation of R(d,m), from the monodromy side).
  Computes, in floating point, with the reduced Burau representation in the convention of Lin-Zhang
  (arXiv:2609.03784, Section 6): for 3 <= d <= 9 and t = zeta_d^r, r = 1, ..., d-1, the invariant
  vector, the signature of the invariant hermitian form on the quotient (r != d/2), the group order
  for d <= 7 by closure of rounded matrices (capped at 300000 elements), the largest eigenvalue
  modulus over random words (fixed seed), and the Galois orbit of r. A sector r counts as finite if
  every Galois conjugate has a definite form (for 10 <= d <= 40 the signatures are taken from
  Lin-Zhang's formula); R(d,m) for 3 <= d <= 12, 1 <= m <= 4 is then counted from the monomials and
  compared with the closed formula. It also lists the d <= 60 for which the Galois orbit of some
  sector lies in {1, d-1}, i.e. for which the conclusion of Lemma 6.1 of Lin-Zhang (stated there for
  d >= 7) fails.
  Expected: "R_dp" = "R_formula" for every (d,m); "lemma61_failures": d = 3, 4, 6 only; group orders
  6 (d = 3), 96 (d = 4, r = 1, 3) and 155520 (d = 6, r = 1, 5); for the other sectors with d <= 7
  the closure exceeds the cap (null).

code/check.py  ->  outputs/check_out.json, outputs/check.out
  Independent checker for producer.py (reads producer_out.json). (a) Exact orders, by closure over
  Z[zeta_3] and Z[i], of the groups generated by the reduced Burau matrices on the quotient, for the
  finite sectors: 6 (d = 3), 96 (d = 4), 155520 (d = 6). (b) The Fermat member, for d <= 9 and
  m in {2, 3}: the number of characters that are Hodge classes by Shioda's criterion (Math. Ann. 245
  (1979)), compared with R(d,m) from producer.py. (c) The signatures of producer.py agree with the
  Chevalley-Weil eigenspace Hodge numbers (r - 1, d - r - 1) of y^d = F(x), up to order, for
  3 <= d <= 9.
  Expected: the orders above; in every "fermat" entry "fermat_hodge_total" >= "R_producer" (as it
  must be, by Corollary 2.4 of the note); "chevalley_weil_signature_match_d3to9": true.

code/burau_sage.py  ->  outputs/burau_sage.json, outputs/burau_sage.out
  Supports: Section 3 (Lemmas 3.1-3.3, Corollary 3.4). Exact arithmetic in Q(zeta_N), N = lcm(4,d).
  For d = 4, 5, 6 and every k = 1, ..., d-1 (q = zeta_d^(-k)), builds McMullen's model
  (Braid groups and Hodge theory, Math. Ann. 355 (2013), Theorem 4.1): the hermitian Gram matrix of
  e_1, ..., e_d (indices mod d) and the maps T_i, and checks: the T_i preserve the form, satisfy the
  braid relations and fix the radical; the radical is spanned by sum e_j and sum q^j e_j; the form on
  the quotient E (dim d - 2) has signature (k - 1, d - 1 - k); T_1, ..., T_(d-1) generate the full
  matrix algebra of E (absolute irreducibility); if some Galois conjugate of k lies in
  [2, d-2], a word w with |sigma(trace w)| > d - 2 for some embedding sigma, so w has infinite order
  (and for q = -1 that T_1 - I is nonzero with square 0); otherwise the order of the group generated
  by the T_i (GAP).
  Sign convention: the script uses <e_i, e_(i+1)> = -sqrt(-1)(1 - q) for the form linear in the
  first variable; see Remark 8.1 of the note for the comparison with the value printed by McMullen.
  Expected: 12 entries, all checks true, "signature_ok": true; group orders 96 (d = 4, k = 1, 3) and
  155520 (d = 6, k = 1, 5); every other entry has an "infinite_order_witness".

code/burau_check.py  ->  outputs/burau_check.json, outputs/burau_check.out
  Supports: Section 3, as a second, independently written implementation: plain Python with its own
  arithmetic in Q(zeta_d), starting from the reduced Burau matrices and the hermitian forms of
  Lin-Zhang, Section 6.3 (for r = d/2 the alternating form times zeta - zeta^(-1)), on the quotient
  by the invariant vector. Same checks as burau_sage.py; group orders by breadth-first closure.
  Expected: 12 entries with all checks true; the same dimensions, irreducibility, infinite-order
  witnesses and group orders as burau_sage.py, and the same signatures up to order (an invariant
  hermitian form on an irreducible representation is unique up to a real factor).

code/compare_outputs.py
  Compares the .json files of a rerun with outputs/ (see Section 2).


4. What the checks do not show
------------------------------
- The checks cover finitely many cases: 3 <= d <= 9 (R(d,m) up to d = 12 in producer.py) and
  1 <= m <= 4 for the dimension counts, d = 4, 5, 6 for the monodromy. The statements for all m are
  proved in the note.
- The monodromy checks concern the two models (McMullen's Theorem 4.1; the reduced Burau
  representation). That these models describe the monodromy of the family of curves y^d = F(u,v) is
  taken from McMullen's paper and Lin-Zhang, as cited in the note.
- An element of infinite order shows that the monodromy group is infinite; that no finite-index
  subgroup has invariant vectors is proved in the note (Lemma 3.3).
- Nothing here concerns the algebraicity of the Hodge classes (Section 6 of the note).
- The group orders of producer.py are computed with rounded floating-point matrices; the exact
  orders are those of check.py and burau_sage.py / burau_check.py.
