# jk.hecke

A simple python package which lets you work with small rank Coxeter
groups and Hecke algebras. Note that we use
[Soergel's normalisation](https://www.ams.org/journals/ert/1997-001-06/S1088-4165-97-00021-6/S1088-4165-97-00021-6.pdf)
where the standard basis satisfies
H<sub>s</sub><sup>2</sup> = 1 + (v<sup>-1</sup>-v)H<sub>s</sup>.

# Usage examples

```python
>>> import jk.hecke.coxeter as c
>>> import jk.hecke.hecke as h
>>> g = c.generate_a3()
>>> h = h.HeckeAlgebra(g)


```