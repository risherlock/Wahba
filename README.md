# Wahba
Implementations of the major developments in 3D-attitude determination of spacecraft.

## Original problem
<p align="center">
  <img src="wahba's_problem.PNG" width="800">
</p>

## Implementations
1. [TRIAD: A Passive System for Determining the Attitude of a Satellite (1964)](matlab/algorithms/triad1964.m)
2. [Davenport's q method: A Vector Approach to the Algebra of Rotations with Applications (1968)](matlab/algorithms/davenport1968.m)
3. [SVD method: Attitude Determination using Vector Observations and Singular Value Decomposition (1968)](matlab/algorithms/svd1968.m)
4. [QUEST: Three-axis Attitude Determination from Vector Observations (1981)](matlab/algorithms/quest1981.m)
5. [FOMA: Attitude Determination using Vector Observations, A Fast Optimal Matrix Algorithm (1993)](matlab/algorithms/foma1993.m)
6. [An Analytic Solution to Wahba's Problem (2013)](matlab/algorithms/yang_analytical2013.m)
7. [Attitude Determination using Newton's Method on Riemannian Manifold (2015)](matlab/algorithms/yang_manifold2015.m)
8. [FLAE: Fast Linear Quaternion Attitude Estimator Using Vector Observations (2017)](matlab/algorithms/flae_newton2017.m) (Newton's method)
9. [FLAE: Fast Linear Quaternion Attitude Estimator Using Vector Observations (2017)](matlab/algorithms/flae_symbolic2017.m) (Symbolic method)
10. [ESOQ, A Closed-Form Solution to the Wahba Problem (1997)](matlab/algorithms/esoq.m)
11. [ESOQ-2 Single-Point Algorithm for Fast Optimal Spacecraft Attitude Determination (1997)](matlab/algorithms/esoq2.m)

## Todos
1. Markley's 12 test conditions
2. Shuster's method of sequential rotation
3. Markley's correction on Shuster's code
4. Handle quaternion double cover
