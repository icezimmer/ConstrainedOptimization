%{
Compute the minimum of a funtion f = x'*Q*x + q'*x in a convex compact domain.
%}

addpath src

% Space dimension and number of simplices
n = 10; K = ceil(0.2*n);
% Kernel dimension and spectral radius of the matrix Q (it must be > 0)
dim_ker = 5; spectral_radius = 10;
% Density of the matrix Q
density = 0.01;
% Minimum value, maximum value and number of zero in the vector q
min_q = 0; max_q = 10; zero_q = 0;
% Seed for the random generator
seed = 4;
% Generate randomly the matrix Q, the vector q and the starting point x_start
[Q, q, P, K_plus, K_avg, date] = GenerateInstance(n, K, dim_ker, spectral_radius, density, min_q, max_q, zero_q, seed);

% Save the parameters
SaveParameters(n, dim_ker, spectral_radius, density, K_plus, K_avg, date)

% Save the variables (not necessary)
SaveVariables(Q, q, P, date)

% List of Frank Wolfe algorithm variants to compare: "Standard", "Exact", "Away-step"
frank_wolfe_variants = ["Exact", "Standard", "Away-step"];

% List of off-the-shelves algorithms by Quadratic Programming to compare with the FW algorithm: "interior-point", "active-set", "sqp"
off_the_shelves = "interior-point";

% Stoping criteria for the algorithms: max error and max number of steps
eps_R = 1e-9; max_steps = 1e3;

% Comparing the methods
[table_results, table_solutions] = ComparingMethods(Q, q, P, date, frank_wolfe_variants, off_the_shelves, eps_R, max_steps);

% Save the results
SaveResults(table_results, table_solutions, date)