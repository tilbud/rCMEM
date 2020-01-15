# MEMs


This repository corresponds to the R package in development for use of Marsh Equilibrium Models.

**Citation:**

_Suggested acknowledgments:_

## Relevant resousrces and literature

See Jim Morris' [MEM web interface for](http://129.252.139.114/model/marsh/mem.asp) inspiration


## Equations and Fundamentals of MEM

Guiding princples:
    - "Equilibrium" in MEM is a target not a description of the current state.  
    - MEM tries to see how simple you can get the model and still be realistic.  
    - MEM ties elevation of marsh surface to ecosystem productivity.  
    - High elevation, stimulates biomass (to a point), and that causes marsh elevation to raise.  
    
### Ideal Mixing Model

A lot of what underpins current versions of the Marsh Equilibrium Model and the Cohort Theory Model is the ideal mixing model, see [@morris2016contributions] and [@holmquist2018accuracy]. This describes the density of a soil as the summation of the soils' relative organic and inorganic fractions, as well as the _self-packing density_ of the organic (k1) and inorganic (k2) matter.

$$BD = {1 \over {{OM\over k1}+ {(1−OM) \over k2}}}$$

If you hold area constant, and you assume you know the input rates of organic and inorganic matter, you can use the ideal mixing model to calculate accretion rate [@morris2016contributions,  @morris2006competition, @morris2007ecological].
