%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain (Minkowski sum of unitary simplices).
%}

addpath src

% Space dimension, number of simplices, force or not non-point-simplices, and the fraction of active constraints respect the solution
n = 1000; K = 10; force_non_point_simplices = true; actv = 0.5;
% Kernel dimension, spectral radius of the matrix Q (considered only if dim_ker<n), and minimum eigenvalue (considered only if dim_ker=0)
dim_ker = 10; spectral_radius = 2; lambda_min = 1;
% Density of the matrix Q
density = 1;

% Seed for the random generator
seed = 2;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, num_vertex, norm_q, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, spectral_radius, lambda_min, density, seed);

% Save the parameters
SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_ker, spectral_radius, lambda_min, density, norm_q, seed, date)

% Save the matrices
SaveMatrices(Q, q, P, date)

% List of Frank Wolfe algorithm variants to compare: "Standard", "Away-step"
frank_wolfe_variants = ["Away-step"];

% List of off-the-shelves algorithms by Quadratic Programming to compare with the FW algorithm: "interior-point-convex", "active-set"
off_the_shelves = ["interior-point-convex"];

% Stoping criteria for the algorithms: max error and max number of steps
eps_RDG = 1e-5; eps_RT = 1e-7; eps_RE = 1e-6; max_steps = 700;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, eps_RDG, eps_RT, eps_RE, max_steps);

% Save the results
SaveBenchmarkResults(table_results, table_solutions, date)