%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain (Minkowski sum of unitary simplices).
%}

addpath src

% Space dimension, number of simplices, force or not non-point-simplices, and the fraction of active constraints respect the solution
n = 100; K = 10; force_non_point_simplices = true; actv = 0.5;
% Kernel dimension, spectral radius of the matrix Q (considered only if dim_ker<n), and minimum eigenvalue (considered only if dim_ker=0)
dim_ker = 10; spectral_radius = 2; lambda_min = 1;
% Density of the matrix Q
density = 1;

% Seed for the random generator
seed = 1;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, num_vertex, norm_q, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, spectral_radius, lambda_min, density, seed);

% Save the parameters
SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_ker, spectral_radius, lambda_min, density, norm_q, date)

% Save the matrices
SaveMatrices(Q, q, P, date)

% List of Frank Wolfe algorithm variants to compare: "Standard", "Away-step"
frank_wolfe_variants = ["Standard", "Away-step"];

% List of off-the-shelves algorithms by Quadratic Programming to compare with the FW algorithm: "interior-point", "active-set", "sqp"
off_the_shelves = [];

% Stoping criteria for the algorithms: max error and max number of steps
eps_R = 1e-10; max_steps = 1e3;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, eps_R, max_steps);

% Save the results
SaveBenchmarkResults(table_results, table_solutions, date)