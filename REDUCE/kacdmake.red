clear kacdet,x;
matrix x(11,11);
array kacdet(6);
for n:=1:6 do begin
for i:=1:11 do for j:=1:11 do x(i,j):=delta(i-j);
for i:=1:p(n,0,0) do for j:=i:p(n,0,0) do      
write x(i,j):=ev((for k:=1:n product l(k)**p(n,i,k))*(for k:=1:n product l(k-n-1)**p(n,j,n+1-k)));
for i:=2:p(n,0,0) do for j:=1:i-1 do x(i,j):=x(j,i);
write kacdet(n):=det(x) 
end;
;end;
