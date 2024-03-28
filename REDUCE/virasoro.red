clear delta,c,h,l,vech;
order h,c,l;
operator delta,l,vech,ev;
noncom l,vech;
for all m let delta(m)=if m=0 then 1 else 0;
for all m,n such that m>n let l(m)*l(n)=l(n)*l(m)+(m-n)*l(m+n)+c*(m**3-m)*delta(m+n)/12;
for all m let vech(m)**2=1;
for all m,n such that m>0 let l(m)*vech(n)=0;
for all m,n such that 0>n let vech(m)*l(n)=0;
for all m let l(0)*vech(m)=h*vech(m),vech(m)*l(0)=vech(m)*h;
for all x let ev(x)=vech(0)*x*vech(0);
;end;
