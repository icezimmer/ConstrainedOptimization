%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and number of simplices
n = 1000; K = 5; force_non_point_simplices = true;
% Kernel dimension and spectral radius of the matrix Q (it must be > 0)
dim_ker = 100; spectral_radius = 10;
% Density of the matrix Q
density = 1;
% Norm of the vector q
% norm_q = 1000;
% Fraction of active constrained for the solution
actv = 0;
% Seed for the random generator
seed = 2;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, dim_ker, force_non_point_simplices, spectral_radius, density, actv, seed);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, date)

% Save the variables (not necessary)
SaveVariables(Q, q, P, date)

% List of Frank Wolfe algorithm variants to compare: "Standard", "Away-step"
frank_wolfe_variants = ["Away-step", "Standard"];

% List of off-the-shelves algorithms by Quadratic Programming to compare with the FW algorithm: "interior-point", "active-set", "sqp"
off_the_shelves = [];%"interior-point";

% Stoping criteria for the algorithms: max error and max number of steps
eps_R = 1e-10; max_steps = 1e4;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, eps_R, max_steps);

% Save the results
SaveResults(table_results, table_solutions, date)