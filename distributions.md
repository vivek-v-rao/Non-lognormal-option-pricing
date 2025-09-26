# Normal distribution

PDF:  
f(x) = (1/(s*sqrt(2*pi))) * exp(-(x-m)^2 / (2*s^2))

Mean = m  
Std dev = s  
Skew = 0  
Kurtosis = 3  

---

# Student t distribution

Parameters: v > 0 (degrees of freedom), location m, scale s > 0

PDF:  
f(x) = [ Gamma((v+1)/2) / (Gamma(v/2)*sqrt(v*pi)*s) ] * (1 + ((x-m)^2)/(v*s^2))^(-(v+1)/2)

Mean = m   (exists if v > 1)  
Std dev = s*sqrt(v/(v-2))   (exists if v > 2)  
Skew = 0   (exists if v > 3)  
Kurtosis = 3*(v-2)/(v-4)   (exists if v > 4)  

---

# Logistic distribution

Parameters: location m, scale s > 0

PDF:  
f(x) = exp(-(x-m)/s) / ( s * (1 + exp(-(x-m)/s))^2 )  
= (1/(4*s)) * sech^2((x-m)/(2*s))

Mean = m  
Std dev = (pi/sqrt(3))*s  
Skew = 0  
Kurtosis = 21/5 = 4.2  

---

# Hyperbolic secant distribution

Parameters: location m, scale s > 0

PDF:  
f(x) = (1/(2*s)) * sech( pi*(x-m)/(2*s) )

Mean = m  
Std dev = s  
Skew = 0  
Kurtosis = 5  

---

# Lognormal distribution

Let ln(X) ~ Normal(m, s^2), with s > 0

PDF (x > 0):  
f(x) = (1/(x*s*sqrt(2*pi))) * exp(-(ln x - m)^2 / (2*s^2))

Mean = exp(m + s^2/2)  
Std dev = sqrt( (exp(s^2)-1) * exp(2*m + s^2) )  
Skew = (exp(s^2)+2) * sqrt(exp(s^2)-1)  
Kurtosis = exp(4*s^2) + 2*exp(3*s^2) + 3*exp(2*s^2) - 3  

---

# Log-logistic distribution

Parameters: scale a > 0, shape b > 0

PDF (x > 0):  
f(x) = (b/a) * (x/a)^(b-1) / (1+(x/a)^b)^2

Mean (exists if b > 1):  
E[X] = a * ( (pi/b) / sin(pi/b) )

Variance (exists if b > 2):  
Var[X] = a^2 * ( (2*pi/b)/sin(2*pi/b) - ((pi/b)/sin(pi/b))^2 )

Skew (exists if b > 3): computed from raw moments  
Kurtosis (exists if b > 4): computed from raw moments  

---

# Log hyperbolic secant distribution

Let ln(X) ~ hyperbolic secant with location m and scale s.  
For |k*s| < pi/2:  
E[X^k] = exp(k*m) * sec(k*s)

PDF (x > 0):  
f(x) = (1/(2*x*s)) * sech( pi*(ln x - m)/(2*s) )

Mean (exists if |s| < pi/2):  
E[X] = exp(m) * sec(s)

Variance (exists if |2*s| < pi/2):  
Var[X] = exp(2*m) * ( sec(2*s) - sec(s)^2 )

Skew (exists if |3*s| < pi/2) and Kurtosis (exists if |4*s| < pi/2) follow from raw moments

