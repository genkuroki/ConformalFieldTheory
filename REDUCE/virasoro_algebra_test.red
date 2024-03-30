centralcharge(tau);
conformalweight(r, s, tau);

partitions(5);
partnum(5);

monomial({1,1,2,3,3,3,4});
monomials(5);
dualmonomial({1,1,2,3,3,3,4});
dualmonomials(5);

lincomb(5);
ltoz(lincomb(5), 5);
allcoeffs(lincomb(5), 5);
solsingvec(2, 2, tau);
singvec(2, 2, tau);

kacmat(3);
kacdet(3);
kacdet_t(3, tau);
kacdet_fact(3, tau);

proj12(singvec(2, 2, tau));
proj12_fact(singvec(2, 2, tau));

act_ff(l(-3), h0, h0+h1-h2);
factorize(act_ff(l(-2)*l(-1), h0, h0+h1-h2));
factorize(act_ff(l(-1)^3, h0, h0+h1-h2));

g_ff(2, 2, cw(1, 1, tau), cw(1, 1, tau) + cw(2, 2, tau) + h2, tau);
f_ff(1, 1, 2, 2, h2, tau);
f_ff_fact(1, 1, 2, 2, h2, tau);
f_ff_sol(1, 1, 2, 2, h2, tau);

;end;
