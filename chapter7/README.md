# Chapter 7

### Jesús Urtasun Elizari - Uiverisity of Milan - 2019/20

Scipy for theoretical physics. Hypergeometrics and polygamma functions.
Reproduce asymptotic expressions in the following paper [https://arxiv.org/pdf/1303.3590.pdf](https://arxiv.org/pdf/1303.3590.pdf)

## Exercise 1 

Import Numpy and Scipy.

```python
import pdb
import numpy as np
from scipy import special as sc
from matplotlib import pyplot as plt
from matplotlib import gridspec
```

Create two files contributions_large_n.py and contributions_small_n.py containing the asymptotic expressions in 

For the large N, implement the expressions in (2.51)

```python
def gamma0(n):
        return (...)

def gamma1(n):
        return (...)

def gamma2(n):
        return (...)
```
For the large N, implement the expressions D (A.6a), Dlog (A.6b), Dhat (A.6c)

```python
# Expressions Dlog (A.6b)
def D0(n):
        return (...)

def D1(n):
        return (...)

def D2(n):
        return (...)

def D3(n):
        return (...)

# Expressions Dlog (A.6b)
def Dlog0(n):
        return - np.log(n)

def Dlog1(n):
        return (1 / 2) * np.log(n) * (2 * np.euler_gamma + np.log(n))

def Dlog2(n):
        return - (1 / 6) * np.log(n) * (6 * pow(np.euler_gamma, 2) + pow(np.pi, 2) + 2 * np.log(n) * (3 * np.euler_gamma + np.log(n)))

def Dlog3(n):
        return sc.polygamma(2, n)

# Expressions Dhat (A.6c)
def Dhat0(n):
        return (...)

def Dhat1(n):
        return (...)

def Dhat2(n):
        return (...)

def Dhat3(n):
        return (...)
```

The special functions such as Rieman zeta and polygamma functions can be called using the scipy.special library
```python
sc.zeta(2)
sc.polygamma(0, n)
```
## Exercise 2 - Small N asymptotics
Implement the small N asymptotics and the subleading contributions
```python
def nlo_small_n_abf(n):

# NLO small N - subleading
def nlo_small_n_abf_subl(n):

# NNLO small N
def nnlo_small_n_abf(n):

# NNLO small N - subleading
def nnlo_small_n_abf_subl(n):

# N3LO small N
def n3lo_small_n_abf(n):

# N3LO small N - subleading
def n3lo_small_n_abf_subl(n):
```

## Exercise 3 - Large N asymptotics

```python
def nlo_large_n(n, exp_type):

        # Soft approximations (2.40)
        csoft = 
        
        return csoft

# NNLO large N
def nnlo_large_n(n, exp_type):

        # Soft approximations (2.40)
        csoft = 

        return csoft

# N3LO large N
def n3lo_large_n(n, exp_type):

        # Soft approximations (2.40)

        csoft = 

        return csoft
```