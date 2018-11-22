# jk.hecke

A simple python package which lets you work with small rank Coxeter
groups and Hecke algebras. Note that we use
[Soergel's normalisation](https://www.ams.org/journals/ert/1997-001-06/S1088-4165-97-00021-6/S1088-4165-97-00021-6.pdf)
where the standard basis satisfies
H<sub>s</sub><sup>2</sup> = 1 + (v<sup>-1</sup>-v)H<sub>s</sup>.

## Notes

The Coxeter identity element is named 'e', and the simple
reflections 'r', 's', 't', etc. We also use the convention that
the identity is the largest element in the left and right orders.

## Usage examples

This example generates the A3 Coxeter group, the corresponding
Hecke algebra, and prints the Kazhdan-Lusztig basis element
corresponding to the identity in the dual Kazhdan-Lusztig basis.

```
>>> import jk.hecke.coxeter as c
>>> import jk.hecke.hecke as h
>>> g = c.generate_a3()
>>> hecke = h.HeckeAlgebra(g)
>>> kl_basis = hecke.generate_kl_basis()
>>> dual_kl_basisa = hecke.generate_dual_kl_basis()
>>> print(kl_basis['e'].dual_kl_filtration())
              e
            r s t
     rs rt sr st ts srts
rsr rst rts srt sts tsr rstsr
   rsrt rsts rtsr srts stsr
      rsrts rstsr srtsr
            rsrtsr
```

This example finds all elements in the same left cell as the
element 'rt'.

```
>>> left_order, right_order = hecke.generate_orders()
Generating digraph
Finished generating digraph
Generating order from digraph
Finished generating order from digraph
>>> for x in g.all_elements():
...     if left_order[x.index][g['rt'].index] and left_order[g['rt'].index][x.index]:
...             print(x)
...
rt
srt
```

More complex example:

```
>>> for x in g.all_elements():
...     if left_order[g['rt'].index][x.inverse().index]:
...             print(f'Multiplying dual KL-basis element rt with KL-basis element {x}:')
...             print((dual_kl_basis['rt'] * kl_basis[x.name]).dual_kl_filtration())
...
Multiplying dual KL-basis element rt with KL-basis element e:
rt

Multiplying dual KL-basis element rt with KL-basis element r:
    rt
t rst rts
    rt

Multiplying dual KL-basis element rt with KL-basis element t:
    rt
r rts tsr
    rt

Multiplying dual KL-basis element rt with KL-basis element rs:
rts
 rt
rts

Multiplying dual KL-basis element rt with KL-basis element rt:
       rt
r t rst rts tsr
  e rs 2rt ts
r t rst rts tsr
       rt

Multiplying dual KL-basis element rt with KL-basis element ts:
rts
 rt
rts

Multiplying dual KL-basis element rt with KL-basis element rst:
    rt
r rts tsr
    rt

Multiplying dual KL-basis element rt with KL-basis element rts:
      rts
    rs rt ts
r t rst 2rts tsr
    rs rt ts
      rts

Multiplying dual KL-basis element rt with KL-basis element tsr:
    rt
t rst rts
    rt
```