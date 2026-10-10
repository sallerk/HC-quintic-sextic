Supporting code and data for the note
"The Hodge conjecture for very general sums of binary quintics and sextics"
==========================================================================

Status: accompanies the note (Kevin Saller, 2026, DOI 10.5281/zenodo.23267886); not peer reviewed.

What the note proves (Theorem A there): for d = 5, 6 and every m >= 1, the Hodge conjecture holds for
X = {F_0(u_0,v_0) + ... + F_m(u_m,v_m) = 0} in P^(2m+1) when F_0, ..., F_m are very general binary forms
of degree d, and dim Hdg^m(X) = 1 + R(d,m) with R(5,m) = 4^(m+1) and
R(6,m) = sum_k (m+1)!/(k!^2 (m+1-2k)!) 16^k 5^(m+1-2k).
The note also treats special members (Theorems B and C: some forms fixed, the others very general
relative to them), integral Hodge classes (Corollary D) and products of independent very general
members (Theorem E).

IMPORTANT: the proofs in the note do not use any of these computations. The scripts below are
consistency checks of the dimension counts (Section 2 of the note), of the monodromy statements
(Section 3) and of the dimension statements for special members (Sections 8 and 9), for finitely
many cases, and of the integral lattice in one case (d = 5, m = 1, Section 10); they are described in
Section 12 ("Computations") of the note. The numbers quoted in the note (values of R(d,m),
h^{2,2} = 1752 for d = 6, m = 2, the dimensions in Remarks 8.6 and 9.6) are evaluations of formulas
stated there; indep_rank.py and fermat_T.py recompute them.

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
burau_check.py 104 s, burau_sage.py 14 s (docker, including start-up). For the scripts added for
Sections 8-10, measured on the same PC, one at a time: dims_gf.py 0.4 s, fermat_T.py 4 s,
i1_dims.py 0.5 s, i3_lines.py 6.5 s, i3_control.py 3.3 s.


2. How to run
-------------
Each script prints a summary to the screen. indep_rank.py, producer.py, check.py, burau_check.py and
burau_sage.py write their result file (.json) into the CURRENT working directory; dims_gf.py,
fermat_T.py, i1_dims.py and i3_lines.py write it to the file named as first argument; i3_control.py
writes no file. The reference results are in outputs/ (the .out files there are the screen output).
To reproduce them, from the top folder of this repository:

    mkdir rerun
    cd rerun
    python ../code/indep_rank.py   > indep_rank.out
    python ../code/producer.py     > producer.out
    python ../code/check.py        > check.out       (reads producer_out.json, so run it in the
                                                      same folder, after producer.py)
    python ../code/burau_check.py  > burau_check.out
    sage -python ../code/burau_sage.py > burau_sage.out
    python ../code/dims_gf.py      dims_gf.json  > dims_gf.out
    python ../code/fermat_T.py     fermat_T.json > fermat_T.out
    python ../code/i1_dims.py      i1_dims.json  > i1_dims.out
    python ../code/i3_lines.py     i3_lines.json > i3_lines.out
    python ../code/i3_control.py   ../code       > i3_control.out   (the argument is the folder
                                                      containing i3_lines.py, which it imports)
    python ../code/compare_outputs.py ../outputs .

With docker instead of a local SageMath, the fifth command is, from the rerun folder (bash; on
Windows with Git Bash, prefix the docker line with MSYS_NO_PATHCONV=1):

    IMG=sagemath/sagemath@sha256:310dfe23fc786f3a5a06a7b0387c60201ec79b32c0db941b7624505ee5edb565
    docker run --rm -v "$(pwd)/..:/w" -w /w/rerun $IMG sage -python ../code/burau_sage.py > burau_sage.out

compare_outputs.py compares every .json file in outputs/ with the file of the same name in the
current folder; it ignores the run times (keys "seconds") and compares floating-point numbers with
relative tolerance 1e-9. Expected last line: "ALL EQUAL". The .out files of burau_sage.py and
burau_check.py contain run times and differ from run to run in those. i3_control.py writes only
i3_control.out; compare it with outputs/i3_control.out directly.

Reproducibility: the files in outputs/ were produced by the commands above in a fresh folder (Windows,
Python and numpy versions as above, SageMath in the docker image above). They agree with the results
of earlier runs of the same code, which has since changed only in comments (compare_outputs.py:
"ALL EQUAL"; four of the five .json files of indep_rank.py, producer.py, check.py, burau_check.py
and burau_sage.py byte for byte, burau_check.json up to run times; the .json files of dims_gf.py,
fermat_T.py, i1_dims.py and i3_lines.py byte for byte). All content is exact integer, rational or
cyclotomic arithmetic or uses fixed random seeds, except the floating-point quantities of producer.py
(residuals "fix_err", "inv_err", eigenvalue moduli "max_abs_eig"), which may differ in the last digits
with other numpy/BLAS versions.


3. Files
--------
Notation as in the note: d = degree of the binary forms, m + 1 = number of blocks, so X has dimension
2m; R = R_0 (x) ... (x) R_m is the Jacobian ring, each factor with Hilbert function
(1 + t + ... + t^(d-2))^2; a multidegree (e_0, ..., e_m) has sector c_s = e_s + 2 (mod d). For a
nonempty T in {0, ..., m}, X_T is the hypersurface defined by the forms F_t, t in T, and T^c is the
complement of T.

code/indep_rank.py  ->  outputs/indep_rank_out.json, outputs/indep_rank.out
  Supports: Section 2 of the note (Lemma 2.2, Definition 2.3 and the closed formula for R(d,m)) and
  the values quoted in Theorem A, Section 12 and Remark 7.4.
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
  Supports: Section 12 (a second computation of R(d,m), from the monodromy side).
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
  braid relations and fix the radical; the radical is spanned by sum e_j and sum q^(-j) e_j
  (McMullen's labelling; the code reverses the cyclic order, see its docstring); the form on
  the quotient E (dim d - 2) has signature (k - 1, d - 1 - k); T_1, ..., T_(d-1) generate the full
  matrix algebra of E (absolute irreducibility); if some Galois conjugate of k lies in
  [2, d-2], a word w with |sigma(trace w)| > d - 2 for some embedding sigma, so w has infinite order
  (and for q = -1 that T_1 - I is nonzero with square 0); otherwise the order of the group generated
  by the T_i (GAP).
  Sign convention: the script uses <e_i, e_(i+1)> = -sqrt(-1)(1 - q) for the form linear in the
  first variable; after reversing the cyclic order of the e_i this is sqrt(-1)(1 - conj(q)), the value
  of McMullen's Theorem 4.1 and of the note.
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

code/dims_gf.py  ->  outputs/dims_gf.json, outputs/dims_gf.out
  Supports: Section 8 (the dimension formula of Theorem 8.4, i.e. of Theorem B) and Section 9
  (Lemma 9.1); described in Section 12, "Special members and products". Written independently of
  indep_rank.py and i1_dims.py.
  Computes the Hilbert series of R with its sector grading as a generating function and checks, for
  3 <= d <= 9 and 1 <= m <= 4: (B) for T = {0, ..., j-1} and every c with c_s = 0 for s not in T,
  the Hodge numbers of V_c(X_F) are (d-1)^|T^c| times those of V_(c_T)(X_T), shifted by
  (|T^c|, |T^c|) (the identity behind the dimension formula of Theorem 8.4); (PSI) for every c with
  all c_s != 0, the Hodge numbers of V_c are those of M_c in Lemma 9.1, i.e. the coefficients of
  prod_s ((<c_s> - 1) + (d - 1 - <c_s>) y), shifted by sigma - 1 with sigma = sum_s <c_s> / d, and
  dim V_c = (d-2)^(m+1).
  Expected: "B_fail": 0, "PSI_fail": 0, "B_checked": 20384, "PSI_checked": 8974. "B_checked"
  counts pairs (T, c): the 17654 trivial ones with T = {0, ..., m} and the 2730 quoted in the note.

code/fermat_T.py  ->  outputs/fermat_T.json, outputs/fermat_T.out
  Supports: Remarks 8.6 and 9.6 of the note; described in Section 12, "Special members and products".
  For F_t = u^d + v^d (t in T) and the other forms very general relative to F_T, computes
  dim Hdg^m(X_F) in two ways. "model": for d = 5 and |T| = j, 1 + 4^(m+1-j) N_j (formula of
  Theorem 8.4), with N_j, the number of primitive Hodge classes of the Fermat quintic of dimension
  2j - 2, from Shioda's criterion; for d = 6 and |T| = 1, 1 + R(6,m) plus the contribution of the
  sectors with c_t in {2, 3, 4}, from the description of Lemma 9.1 with the eigenlines of the Fermat
  sextic curve. "shioda": 1 plus the number of characters of the Fermat hypersurface of degree d and
  dimension 2m that satisfy Shioda's criterion and whose pair sectors c_s = a_(2s) + a_(2s+1) are 0
  for s not in T (d = 5), resp. lie in {0, 1, -1} for s != t (d = 6). Ranges: m = 1, ..., 4, and
  j = 1, ..., m+1 for d = 5.
  Expected: "ok": true in every entry; for d = 5 and j = 2 the values 37, 145, 577, 2305
  (m = 1, ..., 4); for d = 6 the values 58, 702, 8658, 107702, with "very_general" = 1 + R(6,m) =
  58, 606, 6962, 81526 and "all_fermat" (dim Hdg^m of the Fermat sextic) = 86, 1752, 38166 for
  m <= 3 (null for m = 4). For d = 5 the agreement is elementary (a pair of exponents (a, -a) adds d
  to every sum in Shioda's criterion), so there it only checks the programs.

code/i1_dims.py  ->  outputs/i1_dims.json, outputs/i1_dims.out
  Supports: Section 8 (Theorem 8.4, Remark 8.6) and Remark 9.6; described in Section 12, "Special
  members and products".
  (A) The identity (B) of dims_gf.py, sector by sector from the multidegrees of R, for 3 <= d <= 9,
  1 <= m <= 4 and T = {0, ..., k-1}, k = 1, ..., m+1; also that the sector dimensions of X_S add up
  to its primitive Betti number (3 <= d <= 9, 1 <= |S| <= 4). (C) For d = 5, 7 and |T| = 1, 2, 3,
  the dimension 1 + (d-1)^(m+1-|T|) N_prim of Theorem 8.4 for F_T Fermat, with N_prim, the number of
  primitive Hodge classes of the Fermat hypersurface X_T, from Shioda's criterion (for |T| = 1 this
  is 1 + (d-1)^(m+1)). (D) For d = 6 and |T| = 1: the sectors c with c_t in {2, 3, 4} and c_s in
  {0, 1, -1} (s != t) and their Hodge numbers.
  Expected: screen "(A) dimension identity: 28 (d,m) pairs, 0 failures" and "(B) Betti sanity: OK";
  (C) N_prim = 4, 36, 400 (d = 5) and 6, 90, 1860 (d = 7), e.g. "dimHdg" = 37, 145, 577, 2305 for
  d = 5, |T| = 2, m = 1, ..., 4; (D) for m = 2 the two sectors (2, 5, 5) and (4, 1, 1), each with
  (2,2)-part of dimension 48 (Remark 9.6).

code/i3_lines.py  ->  outputs/i3_lines.json, outputs/i3_lines.out
  Supports: Corollary 10.1 (Corollary D) in the case d = 5, m = 1, at the Fermat member; described in
  Section 12, "Special members and products". Plain Python, exact arithmetic.
  With its own arithmetic in Q(zeta_10), computes the intersection matrix of the 75 lines of the
  Fermat quintic surface x_0^5 + x_1^5 + x_2^5 + x_3^5 = 0, the action of G = mu_5^2 (scaling the
  pairs (x_0,x_1) and (x_2,x_3)) on them, the G-invariant part N^G of the lattice N spanned by the
  lines, and the sublattice spanned by the 25 lines Lambda_p (pairing (x_0,x_1) | (x_2,x_3)), with
  their ranks and discriminants. That the 75 lines generate the Neron-Severi group of the Fermat
  quintic surface is taken from Degtyarev (Theorem 1.2 of the paper cited in the note).
  Expected: "rank_all_lines": 37, "rank_split": 17, "rank_G_invariants": 17,
  "disc_NS_fermat": "244140625" (= 5^12), "disc_NS_G_invariant" = "disc_split_span" = "1953125"
  (= 5^9), "index_squared": "1", so the two lattices coincide.

code/i3_control.py  ->  outputs/i3_control.out
  Control for i3_lines.py (imports it from the folder given as argument): doubles one basis vector of
  the lattice spanned by the 25 lines Lambda_p, which must multiply the discriminant by 4 (index 2).
  Expected: "disc split 1953125 disc doubled 7812500 ratio 4 (expect 4)".

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
- Nothing here concerns the algebraicity of the Hodge classes (Sections 6, 8, 9 and 11 of the
  note).
- The group orders of producer.py are computed with rounded floating-point matrices; the exact
  orders are those of check.py and burau_sage.py / burau_check.py.
- dims_gf.py, i1_dims.py and fermat_T.py check dimension identities and counts. fermat_T.py counts
  Hodge classes at the Fermat member; that these counts are the dimensions for the members of
  Remarks 8.6 and 9.6 rests on Proposition 8.3 and Lemma 9.5 of the note, not on a computation.
- i3_lines.py treats only d = 5, m = 1 at the Fermat member, and takes from Degtyarev that the 75
  lines generate the Neron-Severi group; Corollary 10.1 itself rests on Degtyarev-Shimada, as cited in
  the note.
