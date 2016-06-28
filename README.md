# `PPR: Piecewise Polynomial Reconstruction`

<p align="center">
 <img src="../master/img/W-PQM-shear-1.png"> &nbsp &nbsp
 <img src="../master/img/W-PQM-shear-2.png"> &nbsp &nbsp
 <img src="../master/img/W-PQM-shear-3.png">
</p>

The <b>PPR</b> pacakge is a Fortran library for high-order piecewise polynomial reconstruction and conservative integral re-mapping. Such operations can be used to construct high-order finite-volume and/or arbitrary lagrangian-eulerian (ALE) type algorithms for the solution of hyperbolic transport problems.

The <b>PPR</b> package supports a variety of conservative polynomial reconstructions, including the piecewise constant, linear, parabolic and quartic types. Each interpolant can be combined with a selection of slope-limiters, including exact monotonicity-preserving and weighted essential non-oscillatory (WENO) type approaches. Support is provided for both uniform and non-uniform structured grid types. 

`[1]` Darren Engwirda and Maxwell Kelley, A WENO-type slope-limiter for a family of piecewise polynomial methods, http://arxiv.org/abs/1606.08188, 2015. Keywords: Piecewise Parabolic Method (PPM), Piecewise Quartic Method (PQM), Weighted Essentially Non-Oscillatory reconstruction (WENO), Finite-Volume method, Semi-Lagrangian method, Arbitrary Lagrangian-Eulerian method (ALE)
