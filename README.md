## Sage Code Usage

### Computation of the Dimension of the Cuspidal Cohomology

The command `Cuspidal_Cohomology_Dimension(p, Fq)` returns the dimension of the cuspidal cohomology for the level $\Gamma_0(3, p)$ with $\mathbb{F}_q$ coefficients.

INPUT:
* `p` - a rational prime number
* `Fq` - a large finite field of prime order

```python
sage: attach("PATH/cuspidal-cohomology-computations.py")
sage: Cuspidal_Cohomology_Dimension(53, GF(12379))
(53, 2)
```

### Computation of Hecke Operators

The command `Compute_Hecke_Operators(p, Fq, l)` returns the characteristic polynomial of the Hecke operator $E_\ell$.  Additionally, it gives the eigenvalues over $\mathbb{F}_q$ and corresponding eigenvectors.

We note that this characteristic polynomial has $\mathbb{F}_q$ coefficients.  One then must find the appropriate polynomial over $\mathbb{C}$, whose roots are the Hecke eigenvalue $e_\ell$ and its complex conjugate $\bar{e_\ell}$.  

INPUT:
* `p` - a rational prime number
* `Fq` - a large finite field of prime order
* `l` - a rational prime number

```python
sage: attach("PATH/cuspidal-cohomology-computations.py")
sage: Compute_Hecke_Operators(53, GF(12379), 2)
(2,
 T^2 + 4*T + 15,
 [
 (10502, Vector space of degree 2 and dimension 1 over Finite Field of size 12379
 User basis matrix:
 [   1 8856]),
 (1873, Vector space of degree 2 and dimension 1 over Finite Field of size 12379
 User basis matrix:
 [   1 4228])
 ])
```

Note, some additional work is needed following this command 