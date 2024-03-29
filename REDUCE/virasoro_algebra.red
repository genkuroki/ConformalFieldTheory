clear delta, c, h, x, z, l, vec;
order h, c, x, z, l;
operator delta, l, vec;
noncom l, vec;
factor l, vec;
for all m let delta(m) = if m = 0 then 1 else 0;
for all m, n such that m > n let l(m)*l(n) = l(n)*l(m) + (m-n)*l(m+n) + c*(m**3-m)*delta(m+n)/12;
for all h let vec(h)**2 = 1;
for all m, n, h such that m > 0 let l(m)*vec(h) = 0;
for all m, n, h such that 0 > n let vec(h)*l(n) = 0;
for all h let l(0)*vec(h) = h*vec(h), vec(h)*l(0) = vec(h)*h;

procedure centralcharge(tt); 6/tt + 13 + 6*tt;
procedure conformalweight(r, s, tt) = (1-r**2)/(4*tt) + (1-r*s)/2 + (1-s**2)/4*tt;

procedure partitions(n);
  if n = 1 then {{1}} else
    append(
      for each p in partitions(n-1) collect 1 . p,
      for each p in partitions(n-1) join
        if length(p)=1 or first(p) < second(p) then {(first(p)+1) . rest(p)} else {}
    );

procedure monomial(p); for each m in p product l(-m);
procedure monomials(d); for each p in partitions(d) collect monomial(p);
procedure dualmonomial(p); for each m in reverse(p) product l(m);
procedure dualmonomials(d); for each p in partitions(d) collect dualmonomial(p);

procedure lincomb(n); begin
  j := -1;
  return for each ll in monomials(n) sum <<
    j := j + 1;
    if j = 0 then ll else mkid(x, j) * ll
  >>;
end;

procedure ltoz(s, d); begin
  ss := s;
  for m := 1:d do ss := (ss where l(-m) => mkid(z, m));
  return ss;
end;

procedure allcoeffs(s, d); begin
  sz := ltoz(s, d);
  cs := coeff(sz, mkid(z, 1));
  for i := 2:d do
    (cs := for each f in cs join for each g in coeff(f, mkid(z, i)) join if g = 0 then {} else {g});
  return cs;
end;

procedure solsing(r, s, tt); begin
  d := r*s;
  sing := lincomb(d);
  v1 := l(1)*sing*vec(h);
  s1 := ltoz(v1*vec(h), d);
  lineq1 := allcoeffs(s1, d);
  v2 := l(2)*sing*vec(h);
  s2 := ltoz(v2*vec(h), d);
  lineq2 := allcoeffs(s2, d);
  lineq := append(lineq1, lineq2);
  lineq_t := sub(c = centralcharge(tt), h = conformalweight(r, s, tt), lineq);
  xs := for i := 1:(length(partitions(d)) - 1) collect mkid(x, i);
  sol := solve(lineq_t, xs);
  return sol;
end;

procedure sing(r, s, tt); begin
  sing := lincomb(r*s);
  sol := solsing(r, s, tt);
  return sub(first(sol), sing)
end;

on rat;
on div;

;end;
