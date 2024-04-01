on rat;
on div;

% [FF] B. L. Feigin & D. B. Fuchs, Verma modules over the virasoro algebra
% Lecture Notes in Mathematics, Volume 1060, 1984, pp. 230--245,
% Topology,
% General and Algebraic Topology and Applications,
% Proceedings of the International Topological Conference held in Leningrad, August 23-27, 1983

% delta - Knonecker's delta (delta(m) = if m = 0 then 1 else 0)
% c - central charge
% h - conformal weight
% x - unknowns
% z - unknowns
% tau - Use t instead of t because t is the symbol for true in REDUCE
% theta - theta^2 = tau
% l - generators of the Virasoro algebra
% vec - ket |h> or bra <h| vector
% vecproj - L(m)*vevproj() = 0 if m < -2
% gamma - mu + 2h of [FF]
% vecf - w^{-gamma}

clear delta, c, h, x, z, tau, theta, l, vec, vecproj, gamma, vecf;
order h, c, x, z, theta, tau, l;
operator delta, l, vec, vecproj;
noncom l, vec, vec_proj, vecf;
factor l, vec, gamma, vecf;

let theta**2 = tau;
for all m let delta(m) = if m = 0 then 1 else 0;
for all m, n such that m > n let l(m)*l(n) = l(n)*l(m) + (m-n)*l(m+n) + c*(m**3-m)*delta(m+n)/12;
for all h let vec(h)**2 = 1;
for all m, n, h such that m > 0 let l(m)*vec(h) = 0;
for all m, n, h such that 0 > n let vec(h)*l(n) = 0;
for all h let l(0)*vec(h) = h*vec(h), vec(h)*l(0) = vec(h)*h;
let vecproj()**2 = 1;
for all m let vecproj()*l(m) = if m < -2 then 0 else if m = -1 then z1*vecproj() else z2*vecproj();
for all h, gamma, epsilon let vecf(h, gamma)*vecf(h, epsilon) = 1;
for all m, h, gamma let l(m)*vecf(h, gamma) = -(-gamma + h*(m+1)) * vecf(h, gamma - m);

procedure commutator(x, y); x*y - y*x;

procedure centralcharge(tau); 6/tau + 13 + 6*tau;
procedure conformalweight(r, s, tau) = (1-r**2)/(4*tau) + (1-r*s)/2 + (1-s**2)/4*tau;
procedure cw(r, s, tau); conformalweight(r, s, tau);

procedure partitions(n); begin
  if n < 0 then return {};
  if n = 0 then return {{}};
  if n = 1 then return {{1}};
  return append(
    for each p in partitions(n-1) collect 1 . p,
    for each p in partitions(n-1) join
      if length(p)=1 or first(p) < second(p) then {(first(p)+1) . rest(p)} else {}
  );
end;

% partnum(n) = partition number of n
procedure partnum(n); length(partitions(n));

% Define degree by deg(L(-m)) = m
% monomial({a,b,c,...}) = ...*L(-c)*L(-b)*L(-a) (a,b,c,...>0)
% monomials(d) = list of all monomials of degree d
% dualmonomial({a,b,c,...}) = l(a)*L(b)*L(c)*...
% dualmonomials(d) = list of all dualmonomials of degree -d

procedure monomial(p); for each m in p product l(-m);
procedure monomials(d); for each p in partitions(d) collect monomial(p);
procedure dualmonomial(p); for each m in reverse(p) product l(m);
procedure dualmonomials(d); for each p in partitions(d) collect dualmonomial(p);

% kacmat(d) = p(d)xp(d) matrix determinant of which is equal to Kac determinant of degree d
procedure kacmat(d); begin
  scalar basis_r := monomials(d);
  scalar basis_l := dualmonomials(d);
  scalar N := partnum(d);
  scalar j := 0;
  scalar k := 0;
  matrix kac_matrix(N, N);
  for each ll_r in basis_r do <<
    k := k + 1;
    for each ll_l in basis_l do <<
      j := j + 1;
      kac_matrix(j, k) := vec(h) * ll_l * ll_r * vec(h);
    >>;
    j := 0;
  >>;
  return kac_matrix;
end;

% kacdet(d) = Kac determinant of degree d
procedure kacdet(d); begin
  kac_matrix := kacmat(d);
  return det(kac_matrix);
end;

% kacdet_t(d, tau) = Kac determinant of degree d for c = c(tau)
procedure kacdet_t(d, tau); sub(c = centralcharge(tau), kacdet(d));

% kacdet_t_rhs(d, tau) = formula of kacdet_t(d, tau)
procedure kacdet_t_rhs(d, tau); begin
  return for r := 1:d product for s := 1:d product (h - cw(r, s, tau))^partnum(d - r*s);
end;

% kacdet_fact(d, tau) = factorization of kacdet_t(d, tau)
procedure kacdet_fact(d, tau); begin
  scalar kd := kacdet_t(d, tau);
  return {num = factorize(num(kd)), den = den(kd)};
end;

% kacdet_sol(d, tau) = solution of kacdet_t(d, tau) w.r.t. h
procedure kacdet_sol(d, tau); begin
  scalar kd := kacdet_t(d, tau);
  return solve(kd, h)
end;

% lincomb(d) = L(-1)^d + x1 L(-2)*L(-1)^{d-2} + ...
procedure lincomb(d); begin
  scalar j := -1;
  return for each ll in monomials(d) sum <<
    j := j + 1;
    if j = 0 then ll else mkid(x, j) * ll
  >>;
end;

% ltoz(s, d) = replace l(-m) in s by z_m 
procedure ltoz(s, d); begin
  scalar ss := s;
  for m := 1:d do ss := (ss where l(-m) => mkid(z, m));
  return ss;
end;

% allcoeffs(s, d) = all coefficients of s
procedure allcoeffs(s, d); begin
  scalar sz := ltoz(s, d);
  scalar cs := coeff(sz, mkid(z, 1));
  for i := 2:d do
    (cs := for each f in cs join for each g in coeff(f, mkid(z, i)) join if g = 0 then {} else {g});
  return cs;
end;

% solsingvec(r, s, tau) = solution x1, x2, ... of the equation of the singular vector for (r, s)
procedure solsingvec(r, s, tau); begin
  scalar d, sing, v1, s1, lineq1, v2, s2, lineq2, lineq, lineq_t, xs, sol;
  if r = 1 and s = 1 then return {{}};
  d := r*s;
  sing := lincomb(d);
  v1 := l(1)*sing*vec(h);
  s1 := ltoz(v1*vec(h), d);
  lineq1 := allcoeffs(s1, d);
  v2 := l(2)*sing*vec(h);
  s2 := ltoz(v2*vec(h), d);
  lineq2 := allcoeffs(s2, d);
  lineq := append(lineq1, lineq2);
  lineq_t := sub(c = centralcharge(tau), h = conformalweight(r, s, tau), lineq);
  xs := for i := 1:(length(partitions(d)) - 1) collect mkid(x, i);
  sol := solve(lineq_t, xs);
  return sol;
end;

% singvec(r, s, tau) = Virasoro part of the singular vector for (r, s) 
procedure singvec(r, s, tau); begin
  scalar sing := lincomb(r*s);
  scalar sol := solsingvec(r, s, tau);
  return sub(first(sol), sing)
end;

% proj12(s) = replace l(m) (m<-2) in s by 0 and l(-1), l(-2) by z1, z2 respectively
procedure proj12(s); begin
  scalar sp := vecproj() * s;
  return vecproj() * sp;
end;

% proj12_fact(s) = factorization of proj12(s)
procedure proj12_fact(s); begin
  scalar sp := vecproj() * s;
  sp := vecproj() * sp;
  return {num = factorize(num(sp)), den = den(sp)};
end;

% act_ff(s, h, gamma) = result of the action of s on w^{-gamma} with conformal weight h
procedure act_ff(s, h, gamma); begin
  scalar g := s * vecf(h, gamma);
  g := g * vecf(h, gamma);
  return g;
end;

% g_ff(r, s, h, gamma, tau) = result of the action of singvec(r, s, tau) on w^{-gamma} with h
procedure g_ff(r, s, h, gamma, tau); begin
  scalar sing := singvec(r, s, tau);
  return act_ff(sing, h, gamma);
end;

% g_ff_factor(k, l, a, b, h, gamma) = the formula of Theorem 3.2 in [FF]
procedure g_ff_factor(k, l, a, b, h, gamma); begin
  return (gamma^2 
    + ((2*a*(k-a)+k)/theta^2 + (2*b*(l-b)+l)*theta^2 + k*l+k+l-(k-2*a)*(l-2*b)) * gamma
    + ((k-2*a)/theta + (l-2*b)*theta)^2 * h
    + (a/theta + b*theta) * ((a+1)/theta + (b+1)*theta) *
      ((k-a)/theta + (l-b)*theta) * ((k-a+1)/theta + (l-b+1)*theta)
  );
end;

% g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1)
% = another version of the formula in Theorem 3.2 of [FF]
procedure g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1); begin
  scalar k := r1 - 1; 
  scalar l := s1 - 1; 
  scalar a := i1 - 1; 
  scalar b := j1 - 1;
  scalar h := cw(r0, s0, tau);
  scalar gamma := cw(r0, s0, tau) + cw(r1, s1, tau) - h2;
  return g_ff_factor(k, l, a, b, h, gamma);
end;

% g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1)
% = explicit factorization of g_ff_factor_lhs(r0, s0, r1, s1, h2, i1, j1)
procedure g_ff_factor_rhs(r0, s0, r1, s1, h2, i1, j1); begin
  return (
    (h2 - cw(r0+r1+1-2*i1, s0+s1+1-2*j1, tau)) *
    (h2 - cw(r0-(r1+1-2*i1), s0-(s1+1-2*j1), tau))
  );
end;

% f_ff(r1, s1, h0, h1, h2, tau)
% = <h2|Phi(w, singvec_{r1,s1}|h1>)|h0> / w^{h2-h1-h0-r1*s1}
% where c = c(tau)
procedure f_ff(r1, s1, h0, h1, h2, tau); g_ff(r1, s1, h0, h0 + h1 - h2, tau);

% f_fusion(r0, s0, r1, s1, h2, tau)
% = <h2|Phi(w, singvec_{r1,s1}|h1>)|h0> / w^{h2-h1-h0-r1*s1}
% where c = c(tau), h0 = h_{r0,s0}(tau), h1 = h_{r1,s1}(tau). 
procedure f_fusion(r0, s0, r1, s1, h2, tau); begin
  scalar h0 := cw(r0, s0, tau);
  scalar h1 := cw(r1, s1, tau);
  return f_ff(r1, s1, h0, h1, h2, tau);
end;

% f_fusion_rhs(r0, s0, r1, s1, h2, tau)
% = the explicit factorization of f_fusion(r0, s0, r1, s1, h2, tau)
procedure f_fusion_rhs(r0, s0, r1, s1, h2, tau); begin
  return ((-1)^(r1*s1) *
    (for i1 := 1:r1 product for j1 := 1:s1 product h2 - cw(r0+r1+1-2*i1, s0+s1+1-2*j1, tau))
  );
end;

% f_fusion_fact(r0, s0, r1, s1, h2, tau)
% = factorization of f_ff(r1, s1, h0, h1, h2, tau)
procedure f_fusion_fact(r0, s0, r1, s1, h2, tau); begin
  scalar f := f_fusion(r0, s0, r1, s1, h2, tau);
  return {num = factorize(num(f)), den = den(f)};
end;

% f_fusion_fact(r0, s0, r1, s1, h2, tau)
% = solution of f_ff(r1, s1, h0, h1, h2, tau) = 0 w.r.t. h2
procedure f_fusion_sol(r0, s0, r1, s1, h2, tau); begin
  scalar f := f_fusion(r0, s0, r1, s1, h2, tau);
  return solve(f, h2);
end;

% f_minimal_model(p, q) = <h|Phi(w, singvec_{p-1,q-1}|0>)|h> / w^{-r1*s1}
% where c = c(-p/q)
procedure f_minimal_model(p, q); f_ff(p-1, q-1, h, 0, h, -p/q);

% f_minimal_model_fact(p, q) = factorization of f_minimal_model(p, q)
procedure f_minimal_model_fact(p, q); begin
  scalar f := f_minimal_model(p, q);
  scalar topcoeff := first(reverse(coeff(f, h)));
  scalar sol := solve(f, h);
  return {topcoeff, sol};
end;

% f_minimal_model_factor_lhs(p, q)
% = certain specialization of the formula of Theorem 3.2 in [FF].
procedure f_minimal_model_factor_lhs(p, q, r, s, h); begin
  return sub(tau=-p/q, g_ff_factor(p-2, q-2, r-1, s-1, h, 0));
end;

% f_minimal_model_factor_rhs(p, q)
% = explicit simplification of f_minimal_model_factor_lhs(p, q)
procedure f_minimal_model_factor_rhs(p, q, r, w, h); begin
  return -4*(q*r - p*s)^2/(p*q) * (h - cw(r, s, -p/q));
end;

% table_minimal_model(p, q)
% = table of the conformal weights of the minimal model for c=c(-p/q).
procedure table_minimal_model(p, q); begin
  matrix table_mm(p-1, q-1);
  for r := 1:(p-1) do for s := 1:(q-1) do begin
    table_mm(r, s) := cw(r, s, -p/q);
  end;
  return table_mm;
end;

;end;
