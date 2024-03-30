# 2024-03-30 Version

l(m)はVirasoro代数の $L_m$ の交換関係をみなし, 積を取ると自動的にすべての項の積の順序が m が小さい順になるように交換関係が用いられる.

Singular vectors や Kac determinant などの計算が実装されている.

使い方については以下を見て欲しい. 

__警告:__ virasoto_algebta.red の仕様は破壊的に自由に改変されるので, これを読んでいる時点では以下と同様の結果が得られるとは限らないことに注意せよ.

## 実行の記録

Reduce (CSL, rev 6657), 10-Dec-2023 ...

1: in "D:\OneDrive\work\ConformalFieldTheory\REDUCE\virasoro_algebra.red"$

2: centralcharge(tau);

$$
\tau ^{-1}\,\left(6\,\tau ^2+13\,\tau +6\right)
$$

3: conformalweight(r, s, tau);

$$
\tau ^{-1}\,\left(-\frac{1}{4}\,\tau ^2\,s^2+\frac{1}{4}\,\tau ^2-\frac{1}{2}\,\tau \,r\,s+\frac{1}{2}\,\tau -\frac{1}{4}\,r^2+\frac{1}{4}\right)
$$

4: partitions(5);

$$
\left\{\left\{1,1,1,1,1\right\},\left\{1,1,1,2\right\},\left\{1,1,3\right\},\left\{1,2,2\right\},\left\{1,4\right\},\left\{2,3\right\},\left\{5\right\}\right\}
$$

5: partnum(5);

$$
7
$$

6: monomials(5);

$$
\left\{l\left(-1\right)^5,l\left(-2\right)\,l\left(-1\right)^3,l\left(-3\right)\,l\left(-1\right)^2,l\left(-2\right)^2\,l\left(-1\right),l\left(-4\right)\,l\left(-1\right),l\left(-3\right)\,l\left(-2\right),l\left(-5\right)\right\}
$$

7: dualmonomials(5);

$$
\left\{l\left(1\right)^5,l\left(1\right)^3\,l\left(2\right),l\left(1\right)^2\,l\left(3\right),l\left(1\right)\,l\left(2\right)^2,l\left(1\right)\,l\left(4\right),l\left(2\right)\,l\left(3\right),l\left(5\right)\right\}
$$

8: 16: lincomb(3);

$$
l\left(-1\right)^3+\mathit{x}_{1}\,l\left(-2\right)\,l\left(-1\right)+\mathit{x}_{2}\,l\left(-3\right)
$$

9: v1 := l(1)*lincomb(3)*vec(h);

$$
\mathit{v}_{1} := 3\,l\left(-1\right)^2\,\mathit{vec}\left(h\right)\,\left(2\,h+\mathit{x}_{1}+2\right)+2\,l\left(-2\right)\,\mathit{vec}\left(h\right)\,\left(h\,\mathit{x}_{1}+2\,\mathit{x}_{2}\right)
$$

10: s1 := v1*vec(h);

$$
\mathit{s}_{1} := 3\,l\left(-1\right)^2\,\left(2\,h+\mathit{x}_{1}+2\right)+2\,l\left(-2\right)\,\left(h\,\mathit{x}_{1}+2\,\mathit{x}_{2}\right)
$$

11: allcoeffs(s1, 3);

$$
\left\{2\,\left(h\,\mathit{x}_{1}+2\,\mathit{x}_{2}\right),3\,\left(2\,h+\mathit{x}_{1}+2\right)\right\}
$$

12: solve(allcoeffs(s1, 3), {x1, x2});

$$
\left\{\left\{\mathit{x}_{1}=-2\,\left(h+1\right),\mathit{x}_{2}=h\,\left(h+1\right)\right\}\right\}
$$

13: solsingvec(1, 3, tau);

$$
\left\{\left\{\mathit{x}_{1}=4\,\tau ,\mathit{x}_{2}=4\,\tau ^2+2\,\tau \right\}\right\}
$$

14: singvec(1, 3, tau);

$$
l\left(-1\right)^3+4\,\tau \,\left(l\left(-2\right)\,l\left(-1\right)\right)+2\,l\left(-3\right)\,\tau \,\left(2\,\tau +1\right)
$$

15: s13 := singvec(1, 3, tau);

$$
\mathit{s}_{13} := l\left(-1\right)^3+4\,\tau \,\left(l\left(-2\right)\,l\left(-1\right)\right)+2\,l\left(-3\right)\,\tau \,\left(2\,\tau +1\right)
$$

16: l(1) * s13 * vec(conformalweight(1, 3, tau));

$$
0
$$

17: l(2) * s13 * vec(conformalweight(1, 3, tau)) where c => centralcharge(tau);

$$
0
$$

18: kd := kacdet(3);

$$
\mathit{kd} := 48\,h^2\,\left(48\,h^4+22\,h^3\,c-142\,h^3+2\,h^2\,c^2-5\,h^2\,c+102\,h^2+3\,h\,c^2-13\,h\,c-20\,h+c^2+2\,c\right)
$$

19: kd := sub(c=centralcharge(tau), kd);

$$
\begin{aligned}
&
\mathit{kd} := 
144\,\tau ^{-2}\,h^2\times
\\ &
\left(16\,h^4\,\tau ^2+44\,h^3\,\tau ^3+48\,h^3\,\tau ^2+44\,h^3\,\tau +24\,h^2\,\tau ^4+94\,h^2\,\tau ^3+173\,h^2\,\tau ^2+94\,h^2\,\tau +24\,h^2+36\,h\,\tau ^4+\right.
\\ &
\left.130\,h\,\tau ^3+178\,h\,\tau ^2+130\,h\,\tau +36\,h+12\,\tau ^4+56\,\tau ^3+89\,\tau ^2+56\,\tau +12\right)
\end{aligned}
$$

20: den(kd);

$$
\tau ^2
$$

21: factorize(num(kd));

$$
\left\{\left\{144,1\right\},\left\{4\,h\,\tau +2\,\tau +3,1\right\},\left\{h\,\tau +\tau +2,1\right\},\left\{4\,h+3\,\tau +2,1\right\},\left\{h+2\,\tau +1,1\right\},\left\{h,2\right\}\right\}
$$

22: solve(kd, h);

$$
\left\{h=0,h=\tau ^{-1}\,\left(-\tau -2\right),h=-2\,\tau -1,h=\tau ^{-1}\,\left(-\frac{1}{2}\,\tau -\frac{3}{4}\right),h=-\frac{3}{4}\,\tau -\frac{1}{2}\right\}
$$

23: kacdet_t(3, tau);

$$
\begin{aligned}
&
144\,\tau ^{-2}\,h^2\times
\\ &
\left(16\,h^4\,\tau ^2+44\,h^3\,\tau ^3+48\,h^3\,\tau ^2+44\,h^3\,\tau +24\,h^2\,\tau ^4+94\,h^2\,\tau ^3+173\,h^2\,\tau ^2+94\,h^2\,\tau +24\,h^2+36\,h\,\tau ^4+\right.
\\ &
\left.130\,h\,\tau ^3+178\,h\,\tau ^2+130\,h\,\tau +36\,h+12\,\tau ^4+56\,\tau ^3+89\,\tau ^2+56\,\tau +12\right)
\end{aligned}
$$

24: kacdet_fact(3, tau);

$$
\begin{aligned}
\{\mathit{num}=&\{\left\{144,1\right\},\left\{4\,h\,\tau +2\,\tau +3,1\right\},\left\{h\,\tau +\tau +2,1\right\},
\left\{4\,h+3\,\tau +2,1\right\},\left\{h+2\,\tau +1,1\right\},\left\{h,2\right\}\},
\\
\mathit{den}=&\tau ^2\}
\end{aligned}
$$

25: s22 := singvec(2, 2, tau);

$$
\begin{aligned}
&
\mathit{s}_{22} := 
\\ &
l\left(-1\right)^4+\tau ^{-2}\,l\left(-2\right)^2\,\left(\tau ^4-2\,\tau ^2+1\right)+2\,\tau ^{-1}\,l\left(-2\right)\,l\left(-1\right)^2\,\left(\tau ^2+1\right)+2\,\tau ^{-1}\,l\left(-3\right)\,l\left(-1\right)\,\left(\tau ^2+3\,\tau +1\right)+
\\ &
3\,\tau ^{-1}\,l\left(-4\right)\,\left(\tau ^2+2\,\tau +1\right)
\end{aligned}
$$

26: proj12(s22);

$$
\tau ^{-2}\,\left(\tau ^4\,\mathit{z}_{2}^2+2\,\tau ^3\,\mathit{z}_{1}^2\,\mathit{z}_{2}+\tau ^2\,\mathit{z}_{1}^4-2\,\tau ^2\,\mathit{z}_{2}^2+2\,\tau \,\mathit{z}_{1}^2\,\mathit{z}_{2}+\mathit{z}_{2}^2\right)
$$

27: proj12_fact(s22);

$$
\left\{\mathit{num}=\left\{\left\{\tau ^2\,\mathit{z}_{2}+\tau \,\mathit{z}_{1}^2+2\,\tau \,\mathit{z}_{2}+\mathit{z}_{2},1\right\},\left\{\tau ^2\,\mathit{z}_{2}+\tau \,\mathit{z}_{1}^2-2\,\tau \,\mathit{z}_{2}+\mathit{z}_{2},1\right\}\right\},\mathit{den}=\tau ^2\right\}
$$

28: factorize(act_ff(l(-1)^3, h0, h0+h1-h2));

$$
\left\{\left\{\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2}+2,1\right\},\left\{\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2}+1,1\right\},\left\{\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2},1\right\}\right\}
$$

29: factorize(act_ff(l(-2)*l(-1), h0, h0+h1-h2));

$$
\left\{\left\{2\,\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2}+1,1\right\},\left\{\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2},1\right\}\right\}
$$

30: factorize(act_ff(l(-3), h0, h0+h1-h2));

$$
\left\{\left\{3\,\mathit{h}_{0}+\mathit{h}_{1}-\mathit{h}_{2},1\right\}\right\}
$$

31: s12 := singvec(1, 2, tau);

$$
\mathit{s}_{12} := l\left(-1\right)^2+\tau \,l\left(-2\right)
$$

32: factorize(num(sub(c=centralcharge(tau), h0=conformalweight(2, 3, tau), h1=conformalweight(1, 2, tau), act_ff(s12, h0, h0+h1-h2))));

$$
\left\{\left\{15\,\tau ^2+4\,\tau \,\mathit{h}_{2}+14\,\tau +3,1\right\},\left\{3\,\tau ^2+4\,\tau \,\mathit{h}_{2}+6\,\tau +3,1\right\}\right\}
$$

33: solve(num(sub(c=centralcharge(tau), h0=conformalweight(2, 3, tau), h1=conformalweight(1, 2, tau), act_ff(s12, h0, h0+h1-h2))), h2);

$$
\left\{\mathit{h}_{2}=\tau ^{-1}\,\left(-\frac{3}{4}\,\tau ^2-\frac{3}{2}\,\tau -\frac{3}{4}\right),\mathit{h}_{2}=\tau ^{-1}\,\left(-\frac{15}{4}\,\tau ^2-\frac{7}{2}\,\tau -\frac{3}{4}\right)\right\}
$$

34: solve(sub(c=centralcharge(tau), h0=conformalweight(2, 3, tau), h1=conformalweight(1, 2, tau), act_ff(s12, h0, h0+h1-h2)), h2);

$$
\left\{\mathit{h}_{2}=\tau ^{-1}\,\left(-\frac{3}{4}\,\tau ^2-\frac{3}{2}\,\tau -\frac{3}{4}\right),\mathit{h}_{2}=\tau ^{-1}\,\left(-\frac{15}{4}\,\tau ^2-\frac{7}{2}\,\tau -\frac{3}{4}\right)\right\}
$$

35: f_ff_fact(2, 3, 1, 2, h2, tau);

$$
\left\{\mathit{num}=\left\{\left\{4\,\mathit{h}_{2}\,\tau +15\,\tau ^2+14\,\tau +3,1\right\},\left\{4\,\mathit{h}_{2}\,\tau +3\,\tau ^2+6\,\tau +3,1\right\}\right\},\mathit{den}=16\,\tau ^2\right\}
$$

36: f_ff_sol(2, 3, 1, 2, h2, tau);

$$
\left\{\mathit{h}_{2}=\frac{-3\,\tau ^2-6\,\tau -3}{4\,\tau },\mathit{h}_{2}=\frac{-15\,\tau ^2-14\,\tau -3}{4\,\tau }\right\}
$$

37: f_ff_sol(1, 1, 2, 3, h2, tau);

$$
\left\{\mathit{h}_{2}=\frac{2\,\tau +1}{4\,\tau },\mathit{h}_{2}=\frac{-8\,\tau ^2+2\,\tau +1}{4\,\tau },\mathit{h}_{2}=\frac{-8\,\tau ^2-10\,\tau -3}{4\,\tau },\mathit{h}_{2}=\frac{6\,\tau -3}{4\,\tau },\mathit{h}_{2}=\frac{-2\,\tau -3}{4\,\tau }\right\}
$$

38: f_ff_sol(1, 2, 2, 3, h2, tau);

$$
\begin{aligned}
&
\left\{
\mathit{h}_{2}=\frac{\tau ^2+2\,\tau +1}{4\,\tau },
\mathit{h}_{2}=\frac{\tau ^2+2\,\tau -3}{4\,\tau },
\mathit{h}_{2}=\frac{-3\,\tau ^2+2\,\tau +1}{4\,\tau },\right.
\\ &
\quad \left.
\mathit{h}_{2}=\frac{-3\,\tau ^2-6\,\tau -3}{4\,\tau },
\mathit{h}_{2}=\frac{-15\,\tau ^2+2\,\tau +1}{4\,\tau },
\mathit{h}_{2}=\frac{-15\,\tau ^2-14\,\tau -3}{4\,\tau }
\right\}
\end{aligned}
$$




