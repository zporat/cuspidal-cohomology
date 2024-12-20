## Sage Code Usage

The program `cuspidal-cohomology.py` loads two commands; their usages are discussed below.

### Computation of the Dimension of the Cuspidal Cohomology

The command `Cuspidal_Cohomology_Dimension(p, Fq)` returns the dimension of the cuspidal cohomology for the level $\Gamma_0(3, p)$ with $\mathbb{F}_q$ coefficients.

INPUT:
* `p` - a rational prime number
* `Fq` - a large finite field of prime order

```python
sage: attach("PATH/cuspidal-cohomology.py")
sage: Cuspidal_Cohomology_Dimension(53, GF(12379))
(53, 2)
```

### Computation of Hecke Operators

The command `Compute_Hecke_Operators(p, Fq, l_list)` returns the characteristic polynomials of the Hecke operators $E_\ell$ directly on $W/W^\text{nc}$ for each $\ell$ in `l_list`.  Additionally, it gives the eigenvalues over $\mathbb{F}_q$ and corresponding eigenvectors.

We note that this gives the characteristic polynomial of $E_\ell$ on $H^3_{\text{cusp}}(\Gamma, \mathbb{F}_q)$.  One must then lift to the appropriate characteristic polynomial (see section 4.2 of the paper).

INPUT:
* `p` - a rational prime number
* `Fq` - a large finite field of prime order
* `l_list` - a list of rational prime numbers

```python
sage: attach("PATH/cuspidal-cohomology.py")
sage: Compute_Hecke_Operators(53, GF(12379), [2, 3, 5, 7, 11])
2 : T^2 + 4*T + 15 : [15, 4, 1] : 10502 : (1, 8856) : 1873 : (1, 4228)
3 : T^2 + 2*T + 12 : [12, 2, 1] : 10503 : (1, 4228) : 1874 : (1, 8856)
5 : T^2 + 12377*T + 1 : [1, 12377, 1] : 1 : (1, 0) : NONE
7 : T^2 + 6*T + 9 : [9, 6, 1] : 12376 : (1, 0) : NONE
11 : T^2 + 12377*T + 1 : [1, 12377, 1] : 1 : (1, 0) : NONE
```

## Magma Code Usage

The above SageMath code has since been translated into [Magma](http://magma.maths.usyd.edu.au/magma/). The command names and inputs are the same in both programs.  Below are examples of the Magma commands. 

### Computation of the Dimension of the Cuspidal Cohomology
```cpp
> load "PATH/cuspidal-cohomology.magma";
> Cuspidal_Cohomology_Dimension(53, GF(12379));
53 : 2
```

### Computation of Hecke Operators
```cpp
> load "PATH/cuspidal-cohomology.magma";
> Compute_Hecke_Action(53, GF(12379), [2, 3, 5, 7, 11]);
2 : $.1^2 + 4*$.1 + 15 : 10502 : (    1  6666) : 1873 : (    1 11932)
3 : $.1^2 + 2*$.1 + 12 : 1874 : (    1  6666) : 10503 : (    1 11932)
5 : $.1^2 + 12377*$.1 + 1 : 1 : (    1     0) : NONE
7 : $.1^2 + 6*$.1 + 9 : 12376 : (    1     0) : NONE
11 : $.1^2 + 12377*$.1 + 1 : 1 : (    1     0) : NONE
Run complete.
```