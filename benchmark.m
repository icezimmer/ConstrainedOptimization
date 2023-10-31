%{
    Compute the minimum of a quadratic function f = x'*Q*x + q'*x in a convex compact domain 
    (Minkowski sum of unitary simplices)
%}

addpath src

% Force or not non-point-simplices
force_non_point_simplices = true;
% Space dimension, number of simplices, and the fraction of active constraints respect the solution
n = 100; K = 20; actv = 0;
% Kernel dimension, spectral radius of the matrix Q (considered only if dim_ker<n), and minimum positive eigenvalue of Q
dim_ker = 0; spectral_radius = 2; lambda_min = 1;
% Density of the matrix Q
density = 1;

% Seed for the random generator
seed = 1;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, num_vertex, norm_q, date] = GenerateInstance(n, K, force_non_point_simplices, actv, dim_ker, ...
    spectral_radius, lambda_min, density, seed);

% Save the parameters
SaveParameters(n, K_plus, K_avg, num_vertex, actv, dim_ker, spectral_radius, lambda_min, density, norm_q, seed, date)

% Save the matrices
SaveMatrices(Q, q, P, date)

% List of Frank Wolfe algorithm variants to compare (possible choices: "Standard", "Away-step")
frank_wolfe_variants = ["Standard", "Away-step"];

% List of off-the-shelves algorithms by Quadratic Programming to compare with the FW algorithm
% Possible choices: "interior-point-convex", "active-set"
off_the_shelves = [];

% Plot or not the dual comparison of the algorithms
dual_comparison = true;

% Stoping criteria for the algorithms: relative duality gap for the FW,
% relative tollerance for the IPC and max number of steps for both
eps_RDG = 1e-7; eps_RT = 1e-9; max_steps = 1e6;
% Max relative error to consider the solution as correct
eps_RE = 1e-10;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, ...
    dual_comparison, eps_RDG, eps_RT, eps_RE, max_steps);

% Save the results
SaveBenchmarkResults(table_results, table_solutions, date)