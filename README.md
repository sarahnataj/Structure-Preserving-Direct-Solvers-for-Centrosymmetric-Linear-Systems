# Structure-Preserving Direct Solvers for Centrosymmetric Linear Systems

Code for the paper: Numerical Stability of Structure-Preserving Direct Solvers for Centrosymmetric Linear Systems

Authors: Sarah Nataj, Chen Greif, Manfred Trummer

Subjects: direct solution of linear systems, centrosymmetric and skew--centrosymmetric matrices, numerical stability, mixed precision

Pdf link:

Abstract: 	This paper analyzes direct solvers for centrosymmetric linear systems, applying structure-preserving factorizations with a particular focus on assessing their stability. We build on existing algorithms and complement the factorizations with equilibration and mixed-precision computations. The solvers  are applied to linear systems arising from spectral discretizations of partial differential equations and the results demonstrate their effectiveness. We  evaluate the accumulation of roundoff errors during the computational process and their impact on the numerical solution. The study demonstrates that errors originating from the factorization of the matrix and a modified substitution propagate in a stable manner, establishing the direct solver’s robustness. Additionally, we provide insights into the solver's stability by proving a bound on the relative error. 
