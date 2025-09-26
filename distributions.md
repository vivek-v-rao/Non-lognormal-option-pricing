# Normal distribution

**PDF**  
```math
f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
Mean: $\mu$
Std dev: $\sigma$
Skew: 0
Kurtosis: 3

Student t distribution
Parameters: degrees of freedom $\nu > 0$, location $\mu$, scale $\sigma > 0$.

PDF

𝑓
(
𝑥
)
=
Γ
 ⁣
(
𝜈
+
1
2
)
Γ
 ⁣
(
𝜈
2
)
𝜈
𝜋
 
𝜎
(
1
+
(
𝑥
−
𝜇
)
2
𝜈
𝜎
2
)
−
(
𝜈
+
1
)
/
2
f(x)= 
Γ( 
2
ν
​
 ) 
νπ
​
 σ
Γ( 
2
ν+1
​
 )
​
 (1+ 
νσ 
2
 
(x−μ) 
2
 
​
 ) 
−(ν+1)/2
 
Mean: $\mu$ (exists if $\nu > 1$)
Std dev: $\sigma\sqrt{\tfrac{\nu}{\nu-2}}$ (exists if $\nu > 2$)
Skew: 0 (exists if $\nu > 3$)
Kurtosis: $3\frac{\nu-2}{\nu-4}$ (exists if $\nu > 4$)

Logistic distribution
Parameters: location $\mu$, scale $s > 0$.

PDF

𝑓
(
𝑥
)
=
𝑒
−
(
𝑥
−
𝜇
)
/
𝑠
𝑠
(
1
+
𝑒
−
(
𝑥
−
𝜇
)
/
𝑠
)
2
=
1
4
𝑠
 
s
e
c
h
2
 ⁣
(
𝑥
−
𝜇
2
𝑠
)
f(x)= 
s(1+e 
−(x−μ)/s
 ) 
2
 
e 
−(x−μ)/s
 
​
 = 
4s
1
​
 sech 
2
 ( 
2s
x−μ
​
 )
Mean: $\mu$
Std dev: $\tfrac{\pi}{\sqrt{3}}s$
Skew: 0
Kurtosis: $\tfrac{21}{5} = 4.2$

Hyperbolic secant distribution
Parameters: location $\mu$, scale $s > 0$.

PDF

𝑓
(
𝑥
)
=
1
2
𝑠
 
s
e
c
h
 ⁣
(
𝜋
(
𝑥
−
𝜇
)
2
𝑠
)
f(x)= 
2s
1
​
 sech( 
2s
π(x−μ)
​
 )
Mean: $\mu$
Std dev: $s$
Skew: 0
Kurtosis: 5

Lognormal distribution
Let $\ln X \sim \mathcal{N}(\mu,\sigma^2)$, with $\sigma > 0$.

PDF (for $x > 0$)

𝑓
(
𝑥
)
=
1
𝑥
𝜎
2
𝜋
exp
⁡
 ⁣
(
−
(
ln
⁡
𝑥
−
𝜇
)
2
2
𝜎
2
)
f(x)= 
xσ 
2π
​
 
1
​
 exp(− 
2σ 
2
 
(lnx−μ) 
2
 
​
 )
Mean: $\exp(\mu + \tfrac{1}{2}\sigma^2)$
Std dev: $\sqrt{(e^{\sigma^2}-1)e^{2\mu+\sigma^2}}$
Skew: $(e^{\sigma^2}+2)\sqrt{e^{\sigma^2}-1}$
Kurtosis: $e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-3$

Log-logistic distribution
Parameters: scale $\alpha > 0$, shape $\beta > 0$.

PDF (for $x > 0$)

𝑓
(
𝑥
)
=
𝛽
𝛼
(
𝑥
/
𝛼
)
𝛽
−
1
(
1
+
(
𝑥
/
𝛼
)
𝛽
)
2
f(x)= 
α
β
​
  
(1+(x/α) 
β
 ) 
2
 
(x/α) 
β−1
 
​
 
Mean (exists if $\beta > 1$):

𝐸
[
𝑋
]
=
𝛼
𝜋
/
𝛽
sin
⁡
(
𝜋
/
𝛽
)
E[X]=α 
sin(π/β)
π/β
​
 
Variance (exists if $\beta > 2$):

V
a
r
[
𝑋
]
=
𝛼
2
(
2
𝜋
/
𝛽
sin
⁡
(
2
𝜋
/
𝛽
)
−
(
𝜋
/
𝛽
sin
⁡
(
𝜋
/
𝛽
)
)
2
)
Var[X]=α 
2
 ( 
sin(2π/β)
2π/β
​
 −( 
sin(π/β)
π/β
​
 ) 
2
 )
Skew (exists if $\beta > 3$): computed from raw moments.
Kurtosis (exists if $\beta > 4$): computed from raw moments.

Log hyperbolic secant distribution
Let $\ln X \sim \mathrm{sech}(\mu,s)$. Then raw moments exist for $\lvert ks \rvert < \tfrac{\pi}{2}$:

𝐸
[
𝑋
𝑘
]
=
𝑒
𝑘
𝜇
sec
⁡
(
𝑘
𝑠
)
.
E[X 
k
 ]=e 
kμ
 sec(ks).
PDF (for $x > 0$):

𝑓
(
𝑥
)
=
1
2
𝑥
𝑠
 
s
e
c
h
 ⁣
(
𝜋
(
ln
⁡
𝑥
−
𝜇
)
2
𝑠
)
f(x)= 
2xs
1
​
 sech( 
2s
π(lnx−μ)
​
 )
Mean (exists if $\lvert s \rvert < \tfrac{\pi}{2}$):

𝐸
[
𝑋
]
=
𝑒
𝜇
sec
⁡
(
𝑠
)
E[X]=e 
μ
 sec(s)
Variance (exists if $\lvert 2s \rvert < \tfrac{\pi}{2}$):

V
a
r
[
𝑋
]
=
𝑒
2
𝜇
(
sec
⁡
(
2
𝑠
)
−
sec
⁡
2
(
𝑠
)
)
Var[X]=e 
2μ
 (sec(2s)−sec 
2
 (s))
Skew (exists if $\lvert 3s \rvert < \tfrac{\pi}{2}$) and Kurtosis (exists if $\lvert 4s \rvert < \tfrac{\pi}{2}$) follow from the raw-moment formulas.
