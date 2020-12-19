# numerical-analysis
Numerical analysis functions in MATLAB.

### General Usage
From a file in the same directory as NumUtils.m call:
```
NumUtils.MethodName(args)
```

### Methods

 - AB2
	 - Adams-Bashforth 2-step (AB2) LMM for solving IVP.
 - Bisection
	 - Bisection Method.
 - CentralDiff
	 - Central Difference for approximating f'(x0).
 - EstimateFPIter
	 - Estimates the number of iterations required for convergence of Fixed Point Iteration algorithm.
 - EulersMethod
	 - Euler's Method for solving an IVP.
 - FivePointMidpoint
	 - Five-Point Midpoint Formula for approximating f'(x0).
 - FixedPointIter
	 - Fixed Point Iteration.
 - ForwardDiff
	 - Forward Difference for approximating f'(x0).
 - GaussSeidelMethod
	 - Gauss-Seidel method for iteratively solving a linear system of equations.
 - GramSchmidt 
	 - Construct orthogonal polynomials w/ Gram-Schmidt.
 - Jacobian
	 - Symbolically calculates Jacobian for a system.
 - JacobisMethod
	 - Jacobi's method for iteratively solving a linear system of equations
 - Lagrange
	 - Generate Lagrange interpolating polynomial.
 - LLS 
	 - Constructs linear least squares polynomial coefficients
 - LogB
	 - Calculates log(X) with base B.
 - NaturalCubicSpline
	 - Calculates the natural cubic spline for f.
 - NewtonsMethod
	 - Newton's method for root finding problem
 - NewtonsMethodForSystems
	 - Newton's Method for iteratively solving a nonlinear system of equations F(x)=0.
 - QuasiNewton
	 - Quasi-Newton Method for root finding problem using global Bisection Method and local Newton's.
 - QuasiSecant
	 - Quasi-Secant Method for root finding problem using global Bisection Method and local Secant.
 - RK2
	 - Runge-Kutta 2-step (RK2) for solving IVP.
 - SecantMethod
	 - Secant method for root finding problem
 - TaylorPoly
	 - Symbolically calculate the first N terms of a Taylor Polynomial.
 - TaylorPolyNTerm
	 - Symbolically calculate the Nth term of a Taylor Polynomial.
 - ThreePointEndpoint
	 - Three-Point Endpoint Formula for approximating f'(x0).
 - ThreePointMidpoint
	 - Three-Point Midpoint Formula for approximating f'(x0).
 - TruncationError
	 - Symbolically calculate the truncation error associated with Taylor Polynomial approximation.
 - TruncationErrorLagrange
	 - Symbolically calculate the truncation error associated with Lagrange Interpolating Polynomial approximation.

### References

 - Numerical Analysis, 10th Edition by Burden, Faires, and Burden