% L(m)'s satidfy the relation of the Virasoro algebta.

operator O;
noncom O;
commutator(O(m), O(n));

commutator(L(0), L(-2));
commutator(L(0), L(-1));
commutator(L(0), L(0));
commutator(L(0), L(1));
commutator(L(0), L(2));
commutator(L(1), L(-2));
commutator(L(1), L(-1));
commutator(L(1), L(0));
commutator(L(1), L(1));
commutator(L(1), L(2));
commutator(L(2), L(-2));
commutator(L(2), L(-1));
commutator(L(2), L(0));
commutator(L(2), L(1));
commutator(L(2), L(2));

% left action of L(m) on |h>

L(-2)*vec(h);
L(-1)*vec(h);
L(0)*vec(h);
L(1)*vec(h);
L(2)*vec(h);

% right action of L(m) on <h|

vec(h)*L(-2);
vec(h)*L(-1);
vec(h)*L(0);
vec(h)*L(1);
vec(h)*L(2);

% <h|h> = 1

vec(h)*vec(h);

% central charge c(tau)
centralcharge(tau);

% conformal weight h_{r,s}(tau)
conformalweight(r, s, tau);
cw(r, s, tau);

% list of the all partition of n
partitions(5);

% partition number of n
partnum(5);

% monomials and dual monomials of L(m)

monomial({1,1,2,3,3,3,4});
monomials(5);
dualmonomial({1,1,2,3,3,3,4});
dualmonomials(5);

% Kac determinant

kacmat(3);
kacdet(3);
kacdet_t(3, tau);
kacdet_fact(3, tau);
kacdet_sol(3, tau);

kd_rhs := kacdet_t_rhs(3, tau);
kacdet_t(3, tau) / kd_rhs;
{num = factorize(num(kd_rhs)), den = den(kd_rhs)};
solve(kd_rhs, h);

% singular vectors

lincomb(4);
ltoz(lincomb(4), 4);
allcoeffs(lincomb(4), 4);
solsingvec(2, 2, tau);
singvec(2, 2, tau);

% projection of singvec_{r,s}(tau) by L(-1)->z1, L(-2)->z2, L(-3),L(-4)...->0

proj12(singvec(2, 2, tau));
proj12_fact(singvec(2, 2, tau));

proj12_fact(singvec(1, 1, tau));
proj12_fact(singvec(1, 2, tau));
proj12_fact(singvec(1, 3, tau));
proj12_fact(singvec(1, 4, tau));
proj12_fact(singvec(1, 5, tau));
proj12_fact(singvec(1, 6, tau));
proj12_fact(singvec(1, 7, tau));

% action of L(m) on w^{h2-h1-h0}

act_ff(l(-3), h0, h0+h1-h2);
factorize(act_ff(l(-2)*l(-1), h0, h0+h1-h2));
factorize(act_ff(l(-1)^3, h0, h0+h1-h2));

% action of singvec_{r,s} on w^{h2-h1-h0}

g_ff(2, 2, cw(1, 1, tau), cw(1, 1, tau) + cw(2, 2, tau) + h2, tau);

% Verify the factorization of g_ff_factor_lhs

solve(g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1), h2);
solve(g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1), h2);
g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1) / g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1);

% <h2|Phi(w, singvec_{r1,s1}|h1>)|h0> / w^{h2-h1-h0-r1*s1}

f_fusion(1, 1, 2, 2, h2, tau);
f_fusion_fact(1, 1, 2, 2, h2, tau);
f_fusion_sol(1, 1, 2, 2, h2, tau);

f_rhs := f_fusion_rhs(1, 1, 2, 2, h2, tau);
f_fusion(1, 1, 2, 2, h2, tau) / f_rhs;
{num = factorize(num(f_rhs)), den = den(f_rhs)};
solve(f_rhs, h2);

% Verify the explicit factorization of f_fusion(r0, s0, r1, s1, h2, tau)

f_fusion(2, 1, 1, 1, h2, tau) / f_fusion_rhs(2, 1, 1, 1, h2, tau);
f_fusion(2, 2, 1, 2, h2, tau) / f_fusion_rhs(2, 2, 1, 2, h2, tau);
f_fusion(3, 2, 2, 3, h2, tau) / f_fusion_rhs(3, 2, 2, 3, h2, tau);
f_fusion(2, 1, 2, 4, h2, tau) / f_fusion_rhs(2, 1, 2, 4, h2, tau);
f_fusion(1, 1, 3, 3, h2, tau) / f_fusion_rhs(1, 1, 3, 3, h2, tau);
f_fusion(1, 2, 3, 3, h2, tau) / f_fusion_rhs(1, 2, 3, 3, h2, tau);
f_fusion(2, 2, 4, 2, h2, tau) / f_fusion_rhs(2, 2, 4, 2, h2, tau);

% the conformal weights of the minimal model for c=c(-p/q)

f_mm45 := f_minimal_model(4, 5);
factorize(f_mm45);
solve(f_mm45, h);

f_minimal_model_fact(2, 3);
f_minimal_model_fact(2, 5);
f_minimal_model_fact(2, 7);
f_minimal_model_fact(2, 9);
f_minimal_model_fact(2, 11);

f_minimal_model_fact(3, 4);
f_minimal_model_fact(3, 5);
f_minimal_model_fact(3, 7);

f_minimal_model_fact(4, 5);

% Verify explicit simplification of f_minimal_model_factor_lhs(p, q)

f_minimal_model_factor_lhs(p, q, r, s, h);
f_minimal_model_factor_rhs(p, q, r, w, h);
f_minimal_model_factor_lhs(p, q, r, s, h) / f_minimal_model_factor_rhs(p, q, r, s, h);

% table of the conformal weights of the minimal model for c=c(-p/q)

table_minimal_model(2, 3);
table_minimal_model(2, 5);
table_minimal_model(2, 7);
table_minimal_model(2, 9);
table_minimal_model(2, 11);
table_minimal_model(2, 13);

table_minimal_model(3, 4);
table_minimal_model(3, 5);
table_minimal_model(3, 7);
table_minimal_model(3, 8);
table_minimal_model(3, 10);
table_minimal_model(3, 11);
table_minimal_model(3, 13);

table_minimal_model(4, 5);
table_minimal_model(4, 7);
table_minimal_model(4, 9);
table_minimal_model(4, 11);
table_minimal_model(4, 13);

table_minimal_model(5, 6);
table_minimal_model(5, 7);
table_minimal_model(5, 8);
table_minimal_model(5, 9);
table_minimal_model(5, 11);
table_minimal_model(5, 12);
table_minimal_model(5, 13);

;end;
