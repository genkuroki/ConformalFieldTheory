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

L(-2)*vec(h);
L(-1)*vec(h);
L(0)*vec(h);
L(1)*vec(h);
L(2)*vec(h);

vec(h)*L(-2);
vec(h)*L(-1);
vec(h)*L(0);
vec(h)*L(1);
vec(h)*L(2);

vec(h)*vec(h);

centralcharge(tau);
conformalweight(r, s, tau);
cw(r, s, tau);

partitions(5);
partnum(5);

monomial({1,1,2,3,3,3,4});
monomials(5);
dualmonomial({1,1,2,3,3,3,4});
dualmonomials(5);

kacmat(3);
kacdet(3);
kacdet_t(3, tau);
kacdet_fact(3, tau);
kacdet_sol(3, tau);

kd_rhs := kacdet_t_rhs(3, tau);
kacdet_t(3, tau) / kd_rhs;
{num = factorize(num(kd_rhs)), den = den(kd_rhs)};
solve(kd_rhs, h);

lincomb(4);
ltoz(lincomb(4), 4);
allcoeffs(lincomb(4), 4);
solsingvec(2, 2, tau);
singvec(2, 2, tau);

proj12(singvec(2, 2, tau));
proj12_fact(singvec(2, 2, tau));

proj12_fact(singvec(1, 1, tau));
proj12_fact(singvec(1, 2, tau));
proj12_fact(singvec(1, 3, tau));
proj12_fact(singvec(1, 4, tau));
proj12_fact(singvec(1, 5, tau));
proj12_fact(singvec(1, 6, tau));
proj12_fact(singvec(1, 7, tau));

act_ff(l(-3), h0, h0+h1-h2);
factorize(act_ff(l(-2)*l(-1), h0, h0+h1-h2));
factorize(act_ff(l(-1)^3, h0, h0+h1-h2));

g_ff(2, 2, cw(1, 1, tau), cw(1, 1, tau) + cw(2, 2, tau) + h2, tau);

solve(g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1), h2);
solve(g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1), h2);
g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1) / g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1);

f_ff(1, 1, 2, 2, h2, tau);
f_ff_fact(1, 1, 2, 2, h2, tau);
f_ff_sol(1, 1, 2, 2, h2, tau);

f_rhs := f_ff_rhs(1, 1, 2, 2, h2, tau);
f_ff(1, 1, 2, 2, h2, tau) / f_rhs;
{num = factorize(num(f_rhs)), den = den(f_rhs)};
solve(f_rhs, h2);

f_ff(2, 1, 1, 1, h2, tau) / f_ff_rhs(2, 1, 1, 1, h2, tau);
f_ff(2, 2, 1, 2, h2, tau) / f_ff_rhs(2, 2, 1, 2, h2, tau);
f_ff(3, 2, 2, 3, h2, tau) / f_ff_rhs(3, 2, 2, 3, h2, tau);
f_ff(2, 1, 2, 4, h2, tau) / f_ff_rhs(2, 1, 2, 4, h2, tau);
f_ff(1, 1, 3, 3, h2, tau) / f_ff_rhs(1, 1, 3, 3, h2, tau);
f_ff(1, 2, 3, 3, h2, tau) / f_ff_rhs(1, 2, 3, 3, h2, tau);
f_ff(2, 2, 4, 2, h2, tau) / f_ff_rhs(2, 2, 4, 2, h2, tau);

;end;
